/**
 * stats.ts — qPCR ΔΔCt 计算与统计检验（纯 TypeScript，零依赖）
 *
 * 移植自 qPCR-app/app.py（第 81–112 行与第 366–426 行）：
 *   ΔCt = 目标 − 管家；ΔΔCt = ΔCt − 对照组 ΔCt 均值；相对表达 = 2^(−ΔΔCt)
 *   分组均值 ± SD、Welch t 检验、Mann-Whitney U（正态近似+连续性校正）、
 *   BH/Bonferroni 校正、显著性星号规则。
 *
 * 数值算法（Numerical Recipes 标准实现）：
 *   - t 分布 CDF：连分式不完全贝塔函数 betai（对齐 scipy.stats.ttest_ind）
 *   - 正态分布 CDF：正则化不完全伽马函数 gammp（对齐 scipy.stats.norm）
 *
 * 仅使用可擦除 TS 语法（无 enum / namespace / 参数属性），Node 24 可直接运行。
 */

// ---------------------------------------------------------------------------
// 数值基础：不完全贝塔函数（连分式）与不完全伽马函数
// ---------------------------------------------------------------------------

const FPMIN = 1e-300;
const EPS = 3e-14; // 双精度下足以收敛到 ~1e-13 相对精度
const MAXIT = 200;

/** Lanczos 近似的 ln(Γ(x))（Numerical Recipes gammln） */
function gammln(x: number): number {
  const cof = [
    57.1562356658629235, -59.5979603554754912, 14.1360979747417471,
    -0.491913816097620199, 0.339946499848118887e-4, 0.465236289270485756e-4,
    -0.983744753048795646e-4, 0.158088703224912494e-3,
    -0.210264441724104883e-3, 0.217439618115212643e-3,
    -0.16431810653676389e-3, 0.844182239838527433e-4,
    -0.261908384015814087e-4, 0.368991826595316234e-5,
  ];
  let xx = x;
  let tmp = xx + 5.2421875; // 有理逼近的任意常数
  tmp = (xx + 0.5) * Math.log(tmp) - tmp;
  let ser = 0.999999999999997092;
  for (let j = 0; j < 14; j++) {
    xx += 1;
    ser += cof[j] / xx;
  }
  return tmp + Math.log(2.5066282746310005 * ser / x);
}

/** 不完全贝塔函数的连分式求值（Numerical Recipes betacf） */
function betacf(a: number, b: number, x: number): number {
  const qab = a + b;
  const qap = a + 1;
  const qam = a - 1;
  let c = 1;
  let d = 1 - (qab * x) / qap;
  if (Math.abs(d) < FPMIN) d = FPMIN;
  d = 1 / d;
  let h = d;
  for (let m = 1; m <= MAXIT; m++) {
    const m2 = 2 * m;
    let aa = (m * (b - m) * x) / ((qam + m2) * (a + m2));
    d = 1 + aa * d;
    if (Math.abs(d) < FPMIN) d = FPMIN;
    c = 1 + aa / c;
    if (Math.abs(c) < FPMIN) c = FPMIN;
    d = 1 / d;
    h *= d * c;
    aa = (-(a + m) * (qab + m) * x) / ((a + m2) * (qap + m2));
    d = 1 + aa * d;
    if (Math.abs(d) < FPMIN) d = FPMIN;
    c = 1 + aa / c;
    if (Math.abs(c) < FPMIN) c = FPMIN;
    d = 1 / d;
    const del = d * c;
    h *= del;
    if (Math.abs(del - 1.0) < EPS) break;
  }
  return h;
}

/** 正则化不完全贝塔函数 I_x(a, b)（Numerical Recipes betai） */
export function betai(a: number, b: number, x: number): number {
  if (x <= 0) return 0;
  if (x >= 1) return 1;
  const bt = Math.exp(
    gammln(a + b) - gammln(a) - gammln(b) +
    a * Math.log(x) + b * Math.log(1 - x),
  );
  if (x < (a + 1) / (a + b + 2)) {
    return (bt * betacf(a, b, x)) / a;
  }
  return 1 - (bt * betacf(b, a, 1 - x)) / b;
}

/** 级数展开的 P(a, x)（Numerical Recipes gser） */
function gser(a: number, x: number): number {
  const gln = gammln(a);
  if (x <= 0) return 0;
  let ap = a;
  let sum = 1 / a;
  let del = sum;
  for (let n = 0; n < MAXIT; n++) {
    ap += 1;
    del *= x / ap;
    sum += del;
    if (Math.abs(del) < Math.abs(sum) * EPS) break;
  }
  return sum * Math.exp(-x + a * Math.log(x) - gln);
}

/** 连分式的 Q(a, x) = 1 − P(a, x)（Numerical Recipes gcf） */
function gcf(a: number, x: number): number {
  const gln = gammln(a);
  let b = x + 1 - a;
  let c = 1 / FPMIN;
  let d = 1 / b;
  let h = d;
  for (let i = 1; i <= MAXIT; i++) {
    const an = -i * (i - a);
    b += 2;
    d = an * d + b;
    if (Math.abs(d) < FPMIN) d = FPMIN;
    c = b + an / c;
    if (Math.abs(c) < FPMIN) c = FPMIN;
    d = 1 / d;
    const del = d * c;
    h *= del;
    if (Math.abs(del - 1.0) < EPS) break;
  }
  return Math.exp(-x + a * Math.log(x) - gln) * h;
}

/** 正则化下不完全伽马函数 P(a, x)（Numerical Recipes gammp） */
export function gammp(a: number, x: number): number {
  if (x < 0 || a <= 0) return NaN;
  if (x < a + 1) return gser(a, x);
  return 1 - gcf(a, x);
}

/** 标准正态 CDF Φ(z)：Φ(z) = P(1/2, z²/2)（z ≥ 0），利用对称性 */
export function normCdf(z: number): number {
  if (Number.isNaN(z)) return NaN;
  const az = Math.abs(z);
  const p = gammp(0.5, (az * az) / 2);
  return z >= 0 ? 0.5 * (1 + p) : 0.5 * (1 - p);
}

/** 自由度为 df 的 Student t 分布的 P(T > |t|)（双侧尾部之一） */
function tTailAbs(tAbs: number, df: number): number {
  // P(T > t) = 0.5 * I_{df/(df+t²)}(df/2, 1/2)
  return 0.5 * betai(df / 2, 0.5, df / (df + tAbs * tAbs));
}

// ---------------------------------------------------------------------------
// 描述统计
// ---------------------------------------------------------------------------

/** 忽略 NaN 的均值；空数组返回 NaN */
export function mean(xs: number[]): number {
  let s = 0;
  let n = 0;
  for (const x of xs) {
    if (!Number.isNaN(x)) {
      s += x;
      n += 1;
    }
  }
  return n === 0 ? NaN : s / n;
}

/** 忽略 NaN 的样本标准差（ddof=1，与 pandas .std() 一致）；n<2 返回 NaN */
export function sampleSd(xs: number[]): number {
  const vals = xs.filter((x) => !Number.isNaN(x));
  const n = vals.length;
  if (n < 2) return NaN;
  const m = mean(vals);
  let ss = 0;
  for (const x of vals) {
    const d = x - m;
    ss += d * d;
  }
  return Math.sqrt(ss / (n - 1));
}

/** 忽略 NaN 的样本方差（ddof=1） */
function sampleVar(xs: number[]): number {
  const sd = sampleSd(xs);
  return Number.isNaN(sd) ? NaN : sd * sd;
}

// ---------------------------------------------------------------------------
// ΔΔCt 计算
// ---------------------------------------------------------------------------

/** 输入行：至少含 group 列以及 target/hk 两个 Ct 数值列 */
export type QpcrRow = Record<string, unknown>;

/** ΔΔCt 计算结果行（在原行基础上附加 dct / ddct / mrna） */
export interface DdctRow extends QpcrRow {
  dct: number;
  ddct: number;
  mrna: number;
}

export interface DdctResult {
  /** 对照组 ΔCt 均值（忽略 NaN） */
  ncMean: number;
  /** 附加了 dct / ddct / mrna 的行（顺序与输入一致） */
  rows: DdctRow[];
}

/** 宽松数值转换：对齐 pandas.to_numeric(errors="coerce") */
function coerceNum(v: unknown): number {
  if (typeof v === "number") return v;
  if (typeof v === "string") {
    const t = v.trim();
    if (t === "") return NaN;
    const n = Number(t);
    return Number.isNaN(n) ? NaN : n;
  }
  if (v == null) return NaN;
  return NaN;
}

/**
 * 计算 ΔCt / ΔΔCt / 2^(−ΔΔCt)。
 *
 * @param rows      数据行，含 group 列
 * @param targetCol 目标基因 Ct 列名
 * @param hkCol     管家基因 Ct 列名
 * @param ctrlGroup 对照组名（与 group 列字符串比较）
 * @throws 当对照组在 group 列中不存在时抛错（对齐旧版 st.error 行为）
 */
export function computeDdct(
  rows: QpcrRow[],
  targetCol: string,
  hkCol: string,
  ctrlGroup: string,
): DdctResult {
  const groups = rows.map((r) => String(r["group"]));
  if (!groups.includes(ctrlGroup)) {
    throw new Error(`对照组“${ctrlGroup}”未在 group 列中找到。`);
  }
  const dcts = rows.map((r) => coerceNum(r[targetCol]) - coerceNum(r[hkCol]));
  const ctrlDcts = dcts.filter((_, i) => groups[i] === ctrlGroup);
  const ncMean = mean(ctrlDcts);
  const out: DdctRow[] = rows.map((r, i) => {
    const dct = dcts[i];
    const ddct = dct - ncMean;
    return { ...r, dct, ddct, mrna: Math.pow(2, -ddct) };
  });
  return { ncMean, rows: out };
}

// ---------------------------------------------------------------------------
// 分组统计
// ---------------------------------------------------------------------------

export interface GroupStat {
  group: string;
  /** 非 NaN 的 mrna 个数（对齐 pandas agg N="count"） */
  n: number;
  mean: number;
  /** 样本 SD（ddof=1），n<2 时为 NaN */
  sd: number;
}

/**
 * 按 group 分组统计 mrna 的 N / mean / sd。
 * 分组顺序按首次出现顺序（与 pandas groupby(sort=True) 不同；
 * 若需与旧版展示完全一致请按组名排序后使用）。
 */
export function groupStats(rows: DdctRow[]): GroupStat[] {
  const order: string[] = [];
  const buckets = new Map<string, number[]>();
  for (const r of rows) {
    const g = String(r["group"]);
    if (!buckets.has(g)) {
      buckets.set(g, []);
      order.push(g);
    }
    if (!Number.isNaN(r.mrna)) buckets.get(g)!.push(r.mrna);
  }
  return order.map((g) => {
    const vals = buckets.get(g)!;
    return { group: g, n: vals.length, mean: mean(vals), sd: sampleSd(vals) };
  });
}

// ---------------------------------------------------------------------------
// Welch t 检验（双侧，对齐 scipy.stats.ttest_ind(equal_var=False)）
// ---------------------------------------------------------------------------

export interface TTestResult {
  t: number;
  df: number;
  /** 双侧 p 值；任一组 n<2 时返回 NaN */
  p: number;
}

export function welchT(x: number[], y: number[]): TTestResult {
  const xs = x.filter((v) => !Number.isNaN(v));
  const ys = y.filter((v) => !Number.isNaN(v));
  const nx = xs.length;
  const ny = ys.length;
  if (nx < 2 || ny < 2) return { t: NaN, df: NaN, p: NaN };
  const mx = mean(xs);
  const my = mean(ys);
  const vx = sampleVar(xs);
  const vy = sampleVar(ys);
  const sx = vx / nx;
  const sy = vy / ny;
  const se = Math.sqrt(sx + sy);
  if (se === 0) {
    // 两组方差均为 0：均值相等 → p=1，不等 → p=0
    const same = mx === my;
    return { t: same ? 0 : Infinity, df: NaN, p: same ? 1 : 0 };
  }
  const t = (mx - my) / se;
  // Welch–Satterthwaite 自由度
  const df = (sx + sy) * (sx + sy) / ((sx * sx) / (nx - 1) + (sy * sy) / (ny - 1));
  const p = 2 * tTailAbs(Math.abs(t), df);
  return { t, df, p: Math.min(1, p) };
}

// ---------------------------------------------------------------------------
// Mann-Whitney U（正态近似 + 连续性校正 +  ties 校正，
// 对齐 scipy.stats.mannwhitneyu(method="asymptotic", alternative="two-sided")）
// ---------------------------------------------------------------------------

export interface MwuResult {
  /** U1（以第一组计） */
  u: number;
  z: number;
  /** 双侧 p 值；任一组 n<2 时返回 NaN */
  p: number;
}

/** 平均秩（处理 ties），返回与输入等长的秩数组与 ties 修正项 */
function rankData(vals: number[]): { ranks: number[]; tieTerm: number } {
  const idx = vals.map((_, i) => i).sort((a, b) => vals[a] - vals[b]);
  const ranks = new Array<number>(vals.length);
  let tieSum = 0; // Σ(t³ − t)
  let i = 0;
  while (i < vals.length) {
    let j = i;
    while (j + 1 < vals.length && vals[idx[j + 1]] === vals[idx[i]]) j += 1;
    const r = (i + j) / 2 + 1; // 平均秩（1-based）
    for (let k = i; k <= j; k++) ranks[idx[k]] = r;
    const t = j - i + 1;
    if (t > 1) tieSum += t * t * t - t;
    i = j + 1;
  }
  return { ranks, tieTerm: tieSum };
}

export function mannWhitneyU(x: number[], y: number[]): MwuResult {
  const xs = x.filter((v) => !Number.isNaN(v));
  const ys = y.filter((v) => !Number.isNaN(v));
  const n1 = xs.length;
  const n2 = ys.length;
  if (n1 < 2 || n2 < 2) return { u: NaN, z: NaN, p: NaN };
  const all = xs.concat(ys);
  const { ranks, tieTerm } = rankData(all);
  let r1 = 0;
  for (let i = 0; i < n1; i++) r1 += ranks[i];
  const u1 = r1 - (n1 * (n1 + 1)) / 2;
  const n = n1 + n2;
  const mu = (n1 * n2) / 2;
  // scipy: sigma = sqrt(n1*n2/12 * ((n+1) − Σ(t³−t)/(n(n−1))))
  const tieAdj = tieTerm / (n * (n - 1));
  const sigma = Math.sqrt(((n1 * n2) / 12) * (n + 1 - tieAdj));
  if (sigma === 0) return { u: u1, z: NaN, p: NaN };
  const d = u1 - mu;
  // 连续性校正：向均值方向收缩 0.5
  const zc = d === 0 ? 0 : d - Math.sign(d) * 0.5;
  const z = zc / sigma;
  const p = 2 * (1 - normCdf(Math.abs(z)));
  return { u: u1, z, p: Math.min(1, p) };
}

// ---------------------------------------------------------------------------
// 多重检验校正与显著性星号
// ---------------------------------------------------------------------------

/** Benjamini–Hochberg FDR 校正（对齐旧版 bh_fdr；NaN 原样保留） */
export function bhFdr(ps: number[]): number[] {
  const n = ps.length;
  // 对齐 numpy.argsort：NaN 排在最后
  const order = ps.map((_, i) => i).sort((a, b) => {
    const pa = Number.isNaN(ps[a]) ? Infinity : ps[a];
    const pb = Number.isNaN(ps[b]) ? Infinity : ps[b];
    return pa - pb;
  });
  const ranks = new Array<number>(n);
  order.forEach((idx, i) => {
    ranks[idx] = i + 1;
  });
  const adj = ps.map((p, i) => (p * n) / ranks[i]);
  // 反向累计最小值（在排序后的序列上）
  const adjSorted = order.map((idx) => adj[idx]);
  for (let i = n - 2; i >= 0; i--) {
    adjSorted[i] = Math.min(adjSorted[i], adjSorted[i + 1]);
  }
  const out = new Array<number>(n);
  order.forEach((idx, i) => {
    out[idx] = Math.min(adjSorted[i], 1);
  });
  return out;
}

/** Bonferroni 校正：min(p × m, 1) */
export function bonferroni(ps: number[]): number[] {
  const m = ps.length;
  return ps.map((p) => Math.min(p * m, 1));
}

/** 显著性星号规则（对齐旧版 p_to_star） */
export function pToStar(p: number): string {
  if (Number.isNaN(p)) return "n.s.";
  if (p < 1e-4) return "****";
  if (p < 1e-3) return "***";
  if (p < 1e-2) return "**";
  if (p < 0.05) return "*";
  return "n.s.";
}
