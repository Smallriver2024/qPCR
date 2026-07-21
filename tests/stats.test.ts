/**
 * stats.test.ts — 与 scipy 基准值做 1e-10 级对拍 + ΔΔCt 手算验证
 *
 * 基准值由 gen_benchmarks.py（managed python + scipy）生成：
 *   scipy.stats.ttest_ind(equal_var=False)
 *   scipy.stats.mannwhitneyu(alternative="two-sided", method="asymptotic")
 *
 * 运行：node stats.test.ts（Node 24 直接运行可擦除 TS）
 */
import assert from "node:assert/strict";
import {
  computeDdct,
  groupStats,
  welchT,
  mannWhitneyU,
  bhFdr,
  bonferroni,
  pToStar,
} from "../src/lib/stats.ts";

let passed = 0;
function approx(actual: number, expected: number, tol: number, label: string): void {
  const ok =
    (Number.isNaN(actual) && Number.isNaN(expected)) ||
    Math.abs(actual - expected) <= tol;
  assert.ok(
    ok,
    `${label}: actual=${actual} expected=${expected} tol=${tol}`,
  );
  passed += 1;
  console.log(`  ok  ${label} = ${actual}`);
}

// ---------------------------------------------------------------------------
// 1. ΔΔCt 手算验证 + fixtures 端到端
// ---------------------------------------------------------------------------
console.log("[1] computeDdct / ΔΔCt 手算");

const rows = [
  { group: "NC", MYC: 25.1, GAPDH: 20.1 },
  { group: "NC", MYC: 25.3, GAPDH: 20.2 },
  { group: "NC", MYC: 24.9, GAPDH: 20.0 },
  { group: "TreatmentA", MYC: 23.8, GAPDH: 20.1 },
  { group: "TreatmentA", MYC: 23.9, GAPDH: 20.1 },
  { group: "TreatmentA", MYC: 23.7, GAPDH: 20.1 },
  { group: "TreatmentB", MYC: 26.3, GAPDH: 20.1 },
  { group: "TreatmentB", MYC: 26.4, GAPDH: 20.1 },
  { group: "TreatmentB", MYC: 26.2, GAPDH: 20.1 },
];

const { ncMean, rows: out } = computeDdct(rows, "MYC", "GAPDH", "NC");

// NC 三行 ΔCt = 5.0, 5.1, 4.9 → 对照组均值恰为 5.0
approx(ncMean, 5.0, 1e-12, "ncMean (NC 对照组 ΔCt 均值)");
// 第 1 行：ΔCt=5.0, ΔΔCt=0, mrna=1
approx(out[0].dct, 5.0, 1e-12, "NC#1 ΔCt");
approx(out[0].ddct, 0.0, 1e-12, "NC#1 ΔΔCt");
approx(out[0].mrna, 1.0, 1e-12, "NC#1 2^-ΔΔCt");
// 第 4 行（TreatmentA#1）：ΔCt=3.7 → ΔΔCt=-1.3 → 2^1.3 ≈ 2.4623
approx(out[3].dct, 3.7, 1e-12, "TreatmentA#1 ΔCt");
approx(out[3].ddct, -1.3, 1e-12, "TreatmentA#1 ΔΔCt");
approx(out[3].mrna, 2.462288826689834, 1e-12, "TreatmentA#1 2^1.3≈2.4623");
// scipy 基准的完整 mrna 序列
const mrnaExpected = [
  1.0, 0.9330329915368065, 1.0717734625362942,
  2.462288826689834, 2.2973967099940746, 2.6390158215457924,
  0.4352752816480623, 0.40612619817811857, 0.46651649576840437,
];
out.forEach((r, i) =>
  approx(r.mrna, mrnaExpected[i], 1e-12, `mrna[${i}] (scipy 基准)`));

// 对照组不存在 → 抛错
assert.throws(() => computeDdct(rows, "MYC", "GAPDH", "NoSuchGroup"));
passed += 1;
console.log("  ok  对照组缺失时抛错");

// 非数值 Ct 按 NaN 处理（对齐 pandas to_numeric coerce）
const dirty = computeDdct(
  [
    { group: "NC", MYC: "25.1", GAPDH: "20.1" },
    { group: "NC", MYC: "abc", GAPDH: 20.0 },
    { group: "T", MYC: 23.8, GAPDH: 20.1 },
  ],
  "MYC", "GAPDH", "NC",
);
assert.ok(Number.isNaN(dirty.rows[1].mrna));
passed += 1;
console.log("  ok  非数值 Ct → NaN（coerce 行为）");

// ---------------------------------------------------------------------------
// 2. groupStats（scipy/pandas 基准）
// ---------------------------------------------------------------------------
console.log("[2] groupStats");

const gs = groupStats(out);
const gsExpected = [
  { group: "NC", n: 3, mean: 1.0016021513577003, sd: 0.06938411014072068 },
  { group: "TreatmentA", n: 3, mean: 2.466233786076567, sd: 0.17084371914931049 },
  { group: "TreatmentB", n: 3, mean: 0.43597265853152845, sd: 0.03020118808340186 },
];
assert.equal(gs.length, 3);
gsExpected.forEach((e, i) => {
  assert.equal(gs[i].group, e.group);
  assert.equal(gs[i].n, e.n);
  approx(gs[i].mean, e.mean, 1e-12, `${e.group} mean`);
  approx(gs[i].sd, e.sd, 1e-12, `${e.group} sd`);
});

// ---------------------------------------------------------------------------
// 3. Welch t 检验 vs scipy.stats.ttest_ind(equal_var=False)
// ---------------------------------------------------------------------------
console.log("[3] welchT vs scipy");

const vals = {
  NC: [1.0, 0.9330329915368065, 1.0717734625362942],
  TreatmentA: [2.462288826689834, 2.2973967099940746, 2.6390158215457924],
  TreatmentB: [0.4352752816480623, 0.40612619817811857, 0.46651649576840437],
};

const welchBench: Array<[string, string, number, number]> = [
  ["NC", "TreatmentA", -13.75746838043311, 0.0015618066745172565],
  ["NC", "TreatmentB", 12.946634735593861, 0.0015678071300621974],
  ["TreatmentA", "TreatmentB", 20.26896295336467, 0.0018106904264683413],
];
for (const [g1, g2, tExp, pExp] of welchBench) {
  const r = welchT(vals[g1 as keyof typeof vals], vals[g2 as keyof typeof vals]);
  approx(r.t, tExp, 1e-10, `welch t  ${g1} vs ${g2}`);
  approx(r.p, pExp, 1e-10, `welch p  ${g1} vs ${g2}`);
}

// 带 ties 的额外用例（scipy 基准）
const tiesX = [1.0, 2.0, 2.0, 3.5];
const tiesY = [2.0, 3.5, 4.0, 5.0];
const wt = welchT(tiesX, tiesY);
approx(wt.t, -1.8516401995451028, 1e-10, "welch t (ties case)");
approx(wt.p, 0.11530844712391755, 1e-10, "welch p (ties case)");

// 小样本 n<2 → NaN
const small = welchT([1.0], [1.0, 2.0]);
assert.ok(Number.isNaN(small.p));
passed += 1;
console.log("  ok  welch n<2 → NaN");

// ---------------------------------------------------------------------------
// 4. Mann-Whitney U vs scipy (asymptotic, two-sided)
// ---------------------------------------------------------------------------
console.log("[4] mannWhitneyU vs scipy (asymptotic)");

const mwuBench: Array<[string, string, number, number]> = [
  ["NC", "TreatmentA", 0.0, 0.08085559837005224],
  ["NC", "TreatmentB", 9.0, 0.08085559837005224],
  ["TreatmentA", "TreatmentB", 9.0, 0.08085559837005224],
];
for (const [g1, g2, uExp, pExp] of mwuBench) {
  const r = mannWhitneyU(vals[g1 as keyof typeof vals], vals[g2 as keyof typeof vals]);
  approx(r.u, uExp, 1e-12, `mwu U  ${g1} vs ${g2}`);
  approx(r.p, pExp, 1e-10, `mwu p  ${g1} vs ${g2}`);
}

// ties 用例（检验 ties 校正项）
const mt = mannWhitneyU(tiesX, tiesY);
approx(mt.u, 2.5, 1e-12, "mwu U (ties case)");
approx(mt.p, 0.13665824773814744, 1e-10, "mwu p (ties case)");

// 小样本 n<2 → NaN
const msmall = mannWhitneyU([1.0], [1.0, 2.0, 3.0]);
assert.ok(Number.isNaN(msmall.p));
passed += 1;
console.log("  ok  mwu n<2 → NaN");

// ---------------------------------------------------------------------------
// 5. BH / Bonferroni 校正（旧版 bh_fdr / numpy 基准）
// ---------------------------------------------------------------------------
console.log("[5] bhFdr / bonferroni");

const rawPs = welchBench.map(([, , , p]) => p);
const bh = bhFdr(rawPs);
const bf = bonferroni(rawPs);
const bhExpected = [
  0.0018106904264683411, 0.0018106904264683411, 0.0018106904264683411,
];
const bfExpected = [
  0.004685420023551769, 0.004703421390186592, 0.0054320712794050235,
];
bhExpected.forEach((e, i) => approx(bh[i], e, 1e-15, `bh[${i}]`));
bfExpected.forEach((e, i) => approx(bf[i], e, 1e-15, `bonferroni[${i}]`));

// ---------------------------------------------------------------------------
// 6. pToStar 星号规则
// ---------------------------------------------------------------------------
console.log("[6] pToStar");

const starCases: Array<[number, string]> = [
  [NaN, "n.s."],
  [0.00009, "****"],
  [0.0005, "***"],
  [0.005, "**"],
  [0.03, "*"],
  [0.05, "n.s."],
  [0.5, "n.s."],
  [0.0018106904264683411, "**"], // 本 fixtures 的 BH 校正后 p
];
for (const [p, s] of starCases) {
  assert.equal(pToStar(p), s, `pToStar(${p})`);
  passed += 1;
}
console.log("  ok  pToStar 边界全部正确（**** / *** / ** / * / n.s.）");

console.log(`\n全部通过：${passed} 项断言 ✔`);
