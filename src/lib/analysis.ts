/**
 * analysis.ts — ΔΔCt 分析流水线（组装 stats.ts，对齐旧版 app.py 第 304–497 行）
 */
import {
  bhFdr,
  bonferroni,
  computeDdct,
  groupStats,
  mannWhitneyU,
  pToStar,
  welchT,
  type DdctRow,
  type GroupStat,
} from "./stats.ts";

export interface PreparedTable {
  columns: string[];
  rows: Array<Record<string, string | number | null>>;
}

export type TestKind = "welch" | "mwu";
export type PAdjust = "bh" | "bonferroni" | "none";

export const TEST_LABELS: Record<TestKind, string> = {
  welch: "Welch t 检验（近似正态）",
  mwu: "Mann–Whitney U（非参数）",
};

export const ADJUST_LABELS: Record<PAdjust, string> = {
  bh: "Benjamini–Hochberg (FDR)",
  bonferroni: "Bonferroni",
  none: "不校正",
};

export interface AnalysisSettings {
  targetCol: string;
  hkCol: string;
  ctrlGroup: string;
  testKind: TestKind;
  pAdjust: PAdjust;
  /** 显著性横线与柱顶的间距（相对比例），0.02–0.30 */
  yPadding: number;
}

export interface PairRow {
  g1: string;
  g2: string;
  n1: number;
  n2: number;
  test: string;
  pRaw: number;
  pAdj: number;
  signif: string;
}

export interface AnalysisResult {
  /** 原始列 + dct/ddct/mrna */
  columns: string[];
  rows: DdctRow[];
  ncMean: number;
  gstats: GroupStat[];
  /** 按组名排序的分组水平（对齐 pandas groupby(sort=True)） */
  levels: string[];
  means: number[];
  sds: number[];
  groupValues: Record<string, number[]>;
  pairwise: PairRow[];
  settings: AnalysisSettings;
}

// --------------------------------------------------------------------------
// 默认列选择（对齐旧版 pick_default / pick_default_excluding）
// --------------------------------------------------------------------------
function pickDefault(cols: string[], keywords: string[], fallbackIdx = 0): string | null {
  const low = cols.map((c) => c.toLowerCase());
  for (const kw of keywords) {
    for (let i = 0; i < low.length; i++) {
      if (low[i].includes(kw)) return cols[i];
    }
  }
  return cols.length > 0 ? cols[Math.min(fallbackIdx, cols.length - 1)] : null;
}

export function defaultColumns(cols: string[]): { target: string | null; hk: string | null } {
  const hk = pickDefault(cols, ["gapdh", "gap", "actin", "actb", "18s", "reference"], 1);
  const pool = cols.filter((c) => c !== hk);
  const target = pickDefault(pool, ["mc", "target", "gene", "ct"], 0);
  return { target, hk };
}

// --------------------------------------------------------------------------
// 主流水线
// --------------------------------------------------------------------------
export function runAnalysis(table: PreparedTable, s: AnalysisSettings): AnalysisResult {
  if (!table.columns.includes("group")) {
    throw new Error("缺少列：group（区分分组）。请补充后重试。");
  }
  const numCandidates = table.columns.filter((c) => c !== "group" && c !== "sample");
  if (numCandidates.length < 2) {
    throw new Error("检测到的数值列不足（至少需要 2 列 Ct）。");
  }
  if (!s.targetCol || !s.hkCol) throw new Error("请选择目标基因列与管家基因列。");
  if (s.targetCol === s.hkCol) throw new Error("目标基因列与管家基因列不能相同。");

  // group 统一为字符串（对齐旧版 df["group"].astype(str)）
  const rows = table.rows.map((r) => ({ ...r, group: String(r["group"] ?? "") }));
  const { ncMean, rows: work } = computeDdct(rows, s.targetCol, s.hkCol, s.ctrlGroup);

  // 分组统计按组名排序（对齐旧版 groupby 默认 sort=True 的展示）
  const gstats = groupStats(work).sort((a, b) => (a.group < b.group ? -1 : a.group > b.group ? 1 : 0));
  const levels = gstats.map((g) => g.group);
  const means = gstats.map((g) => g.mean);
  const sds = gstats.map((g) => g.sd);

  const groupValues: Record<string, number[]> = {};
  for (const g of levels) groupValues[g] = [];
  for (const r of work) {
    const g = String(r["group"]);
    if (groupValues[g] && !Number.isNaN(r.mrna)) groupValues[g].push(r.mrna);
  }

  // 两两比较
  const pairwise: PairRow[] = [];
  const rawP: number[] = [];
  for (let i = 0; i < levels.length; i++) {
    for (let j = i + 1; j < levels.length; j++) {
      const x = groupValues[levels[i]];
      const y = groupValues[levels[j]];
      let p: number;
      if (x.length < 2 || y.length < 2) {
        p = NaN;
      } else if (s.testKind === "welch") {
        p = welchT(x, y).p;
      } else {
        p = mannWhitneyU(x, y).p;
      }
      rawP.push(p);
      pairwise.push({
        g1: levels[i],
        g2: levels[j],
        n1: x.length,
        n2: y.length,
        test: TEST_LABELS[s.testKind],
        pRaw: p,
        pAdj: NaN, // 校正后回填
        signif: "",
      });
    }
  }
  const adjP =
    s.pAdjust === "bh" ? bhFdr(rawP) : s.pAdjust === "bonferroni" ? bonferroni(rawP) : rawP.slice();
  pairwise.forEach((row, i) => {
    row.pAdj = adjP[i];
    row.signif = pToStar(adjP[i]);
  });

  return {
    columns: [...table.columns, "dct", "ddct", "mrna"],
    rows: work,
    ncMean,
    gstats,
    levels,
    means,
    sds,
    groupValues,
    pairwise,
    settings: s,
  };
}

// --------------------------------------------------------------------------
// 展示格式化
// --------------------------------------------------------------------------
/** 展示保留 4 位小数（对齐旧版 .round(4)，避免浮点尾巴）；非有限值显示空。 */
export function fmt4(v: unknown): string {
  if (v === null || v === undefined) return "";
  if (typeof v === "number") {
    if (!Number.isFinite(v)) return "";
    return String(Math.round(v * 1e4) / 1e4);
  }
  return String(v);
}

/** p 值按 3 位有效数字（对齐 Python f"{p:.3g}"，指数补两位）。 */
export function fmtG3(p: number): string {
  if (Number.isNaN(p)) return "NaN";
  if (p === 0) return "0";
  let s = p.toPrecision(3);
  if (s.includes("e")) {
    const [m, e] = s.split("e");
    const sign = e.startsWith("-") ? "-" : "+";
    const digits = e.replace(/^[+-]/, "").padStart(2, "0");
    s = `${m}e${sign}${digits}`;
  } else if (s.includes(".")) {
    s = s.replace(/\.?0+$/, "");
  }
  return s;
}

/** p 值表格展示：4 位小数，过小显示 <0.0001。 */
export function fmtP(p: number): string {
  if (Number.isNaN(p)) return "";
  if (p < 1e-4) return "<0.0001";
  return String(Math.round(p * 1e4) / 1e4);
}
