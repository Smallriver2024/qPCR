import io, csv
import numpy as np
import pandas as pd
import streamlit as st

st.set_page_config(page_title="qPCR ΔΔCt 计算器", page_icon="🧪", layout="wide")
st.title("qPCR ΔΔCt 在线计算器")
st.caption("复制粘贴或上传数据 → 选择列与对照组 → 一键计算并导出 ExcelR")

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
