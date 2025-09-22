import io, csv
import numpy as np
import pandas as pd
import streamlit as st
import math
import itertools
from typing import List, Tuple
from scipy import stats
import matplotlib.pyplot as plt

st.set_page_config(page_title="qPCR ΔΔCt 计算器", page_icon="🧪", layout="wide")
st.title("qPCR ΔΔCt 在线计算器")
st.caption("复制粘贴或上传数据 → 选择列与对照组 → 一键计算➕绘图➕导出计算结果")

# -------------------- 小工具 --------------------
def read_input(pasted: str, uploaded):
    """优先读上传文件；否则尝试从粘贴文本读。"""
    if uploaded is not None:
        name = uploaded.name.lower()
        if name.endswith((".xlsx", ".xls")):
            return pd.read_excel(uploaded)
        return pd.read_csv(uploaded)
    if not pasted.strip():
        return pd.DataFrame()

    # 尝试自动分隔符（首行探测，失败则回退Tab/逗号）
    try:
        sep = csv.Sniffer().sniff(pasted.splitlines()[0]).delimiter
    except Exception:
        sep = "\t" if "\t" in pasted else ","
    for s in (sep, ",", "\t"):
        try:
            return pd.read_csv(io.StringIO(pasted), sep=s, engine="python")
        except Exception:
            pass
    return pd.DataFrame()

def pick_default(cols, keywords, fallback_idx=0):
    """按关键词优先，找不到用 fallback_idx。"""
    low = [str(c).lower() for c in cols]
    for kw in keywords:
        for i, c in enumerate(low):
            if kw in c:
                return cols[i]
    return cols[min(fallback_idx, len(cols) - 1)] if cols else None

def coerce_num(df, cols):
    for c in cols:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    return df

# -------------------- 侧边栏：数据与参数 --------------------
with st.sidebar:
    st.header("数据与参数")
    mode = st.radio("输入方式", ["复制粘贴", "上传文件"], horizontal=True)
    pasted = st.text_area("粘贴数据（首行为列名）", height=160,
                          placeholder="group\tGAPDH\tGENE\nNC\t12.1\t23.2\nNC\t12.3\t22.8\nTreatmentA\t12.4\t21.9") if mode=="复制粘贴" else ""
    uploaded = st.file_uploader("上传 CSV / Excel", type=["csv","xlsx","xls"]) if mode=="上传文件" else None

    ctrl_group = st.text_input("对照组名称（与 group 列一致）", "NC")
    show_stats = st.checkbox("输出分组均值±SD", True)
    st.divider()
    st.subheader("统计与作图")
    with st.form("stats_form", clear_on_submit=False):
        test_kind = st.selectbox(
            "两两比较检验方式",
            ["Welch t 检验（近似正态）", "Mann–Whitney U（非参数）"],
            index=0
        )
        p_adjust = st.selectbox(
            "多重比较校正",
            ["Benjamini–Hochberg (FDR)", "Bonferroni", "不校正"],
            index=0
        )
        y_padding = st.slider("显著性横线与柱顶的间距（相对比例）", 0.02, 0.30, 0.08, 0.01)
        submitted_stats = st.form_submit_button("更新图表与统计", use_container_width=True)

    run = st.button("计算并生成结果", type="primary", use_container_width=True)


# -------------------- 读取与基础校验 --------------------
df = read_input(pasted, uploaded)
if df.empty:
    st.info("请在左侧粘贴数据或上传文件。至少包含：group 列 + 2 列 Ct 数值（目标与管家）。")
    st.stop()

if "group" not in df.columns:
    st.error("缺少列：group（区分分组）。请补充后重试。")
    st.stop()

st.subheader("数据预览")
st.dataframe(df.head(20), use_container_width=True)

# -------------------- 选择列 --------------------
num_candidates = [c for c in df.columns if c != "group"]
if len(num_candidates) < 2:
    st.error("检测到的数值列不足（至少需要 2 列 Ct）。")
    st.stop()

default_hk = pick_default(num_candidates, ["gapdh","actb","18s","reference"], fallback_idx=1)
default_tg = pick_default(num_candidates, ["mc","target","gene","ct"], fallback_idx=0)

col1, col2 = st.columns(2)
with col1:
    target_col = st.selectbox("目标基因 Ct 列", options=num_candidates,
                              index=num_candidates.index(default_tg) if default_tg in num_candidates else 0)
with col2:
    hk_col = st.selectbox("管家基因 Ct 列", options=num_candidates,
                          index=num_candidates.index(default_hk) if default_hk in num_candidates else min(1, len(num_candidates)-1))

# -------------------- 计算与导出 --------------------
if not run:
    st.info("设置好对照组与列名后，点击左侧按钮开始计算。")
    st.stop()

work = df.copy()
work = coerce_num(work, [target_col, hk_col])
work["dct"] = work[target_col] - work[hk_col]

groups = work["group"].astype(str)
if ctrl_group not in groups.unique():
    st.error(f"对照组“{ctrl_group}”未在 group 列中找到。")
    st.stop()

nc_mean = work.loc[groups == ctrl_group, "dct"].mean(skipna=True)
work["ddct"] = work["dct"] - nc_mean
work["mrna"] = np.power(2.0, -work["ddct"])

st.success("计算完成")
st.dataframe(work, use_container_width=True)

group_stats = None
if show_stats:
    group_stats = (work.groupby("group", dropna=False)["mrna"]
                        .agg(N="count", mrna_mean="mean", mrna_sd="std")
                        .reset_index())
    st.markdown("**分组均值 ± SD（mrna）**")
    st.dataframe(group_stats, use_container_width=True)

buf = io.BytesIO()
with pd.ExcelWriter(buf, engine="openpyxl") as w:
    work.to_excel(w, index=False, sheet_name="qpcr_results")
    pd.DataFrame({
        "Parameter": ["Control group", "Target Ct column", "Housekeeper Ct column", "nc_mean(ΔCt)"],
        "Value":     [ctrl_group, target_col, hk_col, nc_mean]
    }).to_excel(w, index=False, sheet_name="settings")
    if show_stats and group_stats is not None:
        group_stats.to_excel(w, index=False, sheet_name="group_stats")

st.download_button("⬇️ 下载 Excel 结果", buf.getvalue(),
                   file_name="qpcr_ddct_results.xlsx",
                   mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                   use_container_width=True)

# -------------------- 柱状图 + 统计检验（自动识别分组） --------------------
# 准备数据：按分组自动聚合
levels: List[str] = list(group_stats["group"].astype(str))
means = group_stats["mrna_mean"].to_numpy()
sds   = group_stats["mrna_sd"].to_numpy()
ns    = group_stats["N"].astype(int).to_numpy()

# 分组对应的观测值列表，用于两两统计
group_values = {g: work.loc[work["group"].astype(str) == g, "mrna"].dropna().to_numpy()
                for g in levels}

# ---- 两两比较：原始 p 值 ----
pairs: List[Tuple[str, str]] = []
raw_pvals: List[float] = []

for (g1, g2) in itertools.combinations(levels, 2):
    x = group_values[g1]
    y = group_values[g2]
    if len(x) < 2 or len(y) < 2:
        # 样本数太小，跳过或给 NaN
        p = np.nan
    else:
        if test_kind.startswith("Welch"):
            # Welch t-test（不等方差更稳健）
            _, p = stats.ttest_ind(x, y, equal_var=False, nan_policy="omit")
        else:
            # Mann-Whitney U（双尾）
            try:
                _, p = stats.mannwhitneyu(x, y, alternative="two-sided")
            except ValueError:
                p = np.nan
    pairs.append((g1, g2))
    raw_pvals.append(p)

raw_pvals = np.array(raw_pvals, dtype=float)

# ---- 多重比较校正 ----
def bh_fdr(p):
    """Benjamini–Hochberg FDR 校正（不依赖 statsmodels）"""
    p = np.asarray(p, dtype=float)
    n = p.size
    order = np.argsort(p)
    ranks = np.empty_like(order)
    ranks[order] = np.arange(1, n + 1)
    adj = p * n / ranks
    # 保证单调不增
    adj_sorted = np.minimum.accumulate(adj[order][::-1])[::-1]
    adj_final = np.empty_like(adj_sorted)
    adj_final[order] = np.minimum(adj_sorted, 1.0)
    return adj_final

if p_adjust.startswith("Benjamini"):
    adj_pvals = bh_fdr(raw_pvals)
elif p_adjust.startswith("Bonferroni"):
    adj_pvals = np.minimum(raw_pvals * len(raw_pvals), 1.0)
else:
    adj_pvals = raw_pvals.copy()

# ---- 显著性标记文本 ----
def p_to_star(p):
    if np.isnan(p):
        return "n.s."
    if p < 1e-4:
        return "****"
    if p < 1e-3:
        return "***"
    if p < 1e-2:
        return "**"
    if p < 0.05:
        return "*"
    return "n.s."

# ---- 画图（均值±SD 柱状图 + 显著性横线） ----
fig, ax = plt.subplots(figsize=(5, 5), dpi=300)
xpos = np.arange(len(levels))

# 柱状图与误差线
bar_container = ax.bar(xpos, means, yerr=sds, capsize=5, alpha=0.9, linewidth=0.6, edgecolor="black")
ax.set_xticks(xpos, levels, rotation=0)
ax.set_ylabel("Relative mRNA Expression")
ax.set_title("qPCR group comparison")

# y 轴上限，留出给显著性横线
ymax_data = np.nanmax(means + sds) if np.isfinite(np.nanmax(means + sds)) else np.nanmax(means)
if not np.isfinite(ymax_data):
    ymax_data = 1.0
ax.set_ylim(0, ymax_data * (1.0 + max(0.25, y_padding * (len(pairs)))))

# 为了避免横线重叠：用一个“占位层”算法安放高度
used_spans: List[Tuple[float, float, float]] = []  # (x_left, x_right, y_level)
def find_free_y(xl, xr, base):
    """给区间 [xl, xr] 找一个不冲突的y层"""
    h = base
    bump = ymax_data * y_padding
    while any((xl <= uxr and xr >= uxl and abs(h - uy) < bump*0.9) for (uxl, uxr, uy) in used_spans):
        h += bump
    used_spans.append((xl, xr, h))
    return h

# 逐对添加显著性横线与标注
for (i, (g1, g2)) in enumerate(pairs):
    j1, j2 = levels.index(g1), levels.index(g2)
    x1, x2 = min(j1, j2), max(j1, j2)

    # 该对比较的基线：两个柱顶的较高者
    top_y = max(means[x1] + (sds[x1] if np.isfinite(sds[x1]) else 0),
                means[x2] + (sds[x2] if np.isfinite(sds[x2]) else 0))
    base_y = top_y + ymax_data * y_padding * 0.4

    y = find_free_y(x1 - 0.2, x2 + 0.2, base_y)

    # 画横线
    ax.plot([x1, x1, x2, x2], [y, y + ymax_data * 0.01, y + ymax_data * 0.01, y],
            linewidth=1.0)

    # 标注文本（采用校正后的 p）
    label = p_to_star(adj_pvals[i]) + ("" if np.isnan(raw_pvals[i]) else f"\n(p={adj_pvals[i]:.3g})")
    ax.text((x1 + x2) / 2, y + ymax_data * 0.012, label,
            ha="center", va="bottom", fontsize=9)

st.pyplot(fig, use_container_width=True)

# ---- 导出两两比较结果表 ----
comp_df = pd.DataFrame({
    "Group1": [a for (a, _) in pairs],
    "Group2": [b for (_, b) in pairs],
    "N1":     [len(group_values[a]) for (a, _) in pairs],
    "N2":     [len(group_values[b]) for (_, b) in pairs],
    "Test":   [test_kind] * len(pairs),
    "p_raw":  raw_pvals,
    "p_adj":  adj_pvals,
    "Signif": [p_to_star(p) for p in adj_pvals],
})
st.markdown("**两两比较结果**")
st.dataframe(comp_df, use_container_width=True)

# 附加到 Excel 导出
buf_sig = io.BytesIO()
with pd.ExcelWriter(buf, engine="openpyxl", mode="a", if_sheet_exists="replace") as w:
    comp_df.to_excel(w, index=False, sheet_name="pairwise_tests")

# 也提供单独下载 CSV（可选）
csv_bytes = comp_df.to_csv(index=False).encode("utf-8")
st.download_button("⬇️ 下载两两比较CSV", csv_bytes, file_name="qpcr_pairwise_tests.csv",
                   mime="text/csv", use_container_width=True)
