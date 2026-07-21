/**
 * exporters.ts — 导出：Excel（4 sheet）/ 两两比较 CSV（BOM）/ 柱状图 PNG（2x）
 */
import * as XLSX from "xlsx";
import { ADJUST_LABELS, fmtG3, TEST_LABELS, type AnalysisResult } from "./analysis.ts";

function downloadBlob(filename: string, blob: Blob): void {
  const url = URL.createObjectURL(blob);
  const a = document.createElement("a");
  a.href = url;
  a.download = filename;
  document.body.appendChild(a);
  a.click();
  document.body.removeChild(a);
  setTimeout(() => URL.revokeObjectURL(url), 1000);
}

/** 单元格值：保留全精度数字，null/NaN 导出为空。 */
function cell(v: unknown): string | number {
  if (v === null || v === undefined) return "";
  if (typeof v === "number") return Number.isFinite(v) ? v : "";
  return String(v);
}

/** Excel：qpcr_results / settings / group_stats / pairwise_tests 四个 sheet。 */
export function exportExcel(result: AnalysisResult): void {
  const s = result.settings;
  const wb = XLSX.utils.book_new();

  const resultsAoa = [
    result.columns,
    ...result.rows.map((r) => result.columns.map((c) => cell(r[c]))),
  ];
  XLSX.utils.book_append_sheet(wb, XLSX.utils.aoa_to_sheet(resultsAoa), "qpcr_results");

  const settingsAoa = [
    ["Parameter", "Value"],
    ["Control group", s.ctrlGroup],
    ["Target Ct column", s.targetCol],
    ["Housekeeper Ct column", s.hkCol],
    ["nc_mean(ΔCt)", cell(result.ncMean)],
    ["Test", TEST_LABELS[s.testKind]],
    ["P adjustment", ADJUST_LABELS[s.pAdjust]],
  ];
  XLSX.utils.book_append_sheet(wb, XLSX.utils.aoa_to_sheet(settingsAoa), "settings");

  const gstatsAoa = [
    ["group", "N", "mrna_mean", "mrna_sd"],
    ...result.gstats.map((g) => [g.group, g.n, cell(g.mean), cell(g.sd)]),
  ];
  XLSX.utils.book_append_sheet(wb, XLSX.utils.aoa_to_sheet(gstatsAoa), "group_stats");

  const pairwiseAoa = [
    ["Group1", "Group2", "N1", "N2", "Test", "p_raw", "p_adj", "Signif"],
    ...result.pairwise.map((p) => [
      p.g1,
      p.g2,
      p.n1,
      p.n2,
      p.test,
      cell(p.pRaw),
      cell(p.pAdj),
      p.signif,
    ]),
  ];
  XLSX.utils.book_append_sheet(wb, XLSX.utils.aoa_to_sheet(pairwiseAoa), "pairwise_tests");

  const out = XLSX.write(wb, { type: "array", bookType: "xlsx" }) as ArrayBuffer;
  downloadBlob(
    "qpcr_ddct_results.xlsx",
    new Blob([out], {
      type: "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
    }),
  );
}

/** 两两比较 CSV（utf-8-sig，对齐旧版 to_csv + BOM）。 */
export function exportCsv(result: AnalysisResult): void {
  const esc = (v: string): string =>
    /[",\n]/.test(v) ? `"${v.replace(/"/g, '""')}"` : v;
  const num = (p: number): string => (Number.isNaN(p) ? "" : fmtG3(p));
  const lines = [
    "Group1,Group2,N1,N2,Test,p_raw,p_adj,Signif",
    ...result.pairwise.map((p) =>
      [esc(p.g1), esc(p.g2), p.n1, p.n2, esc(p.test), num(p.pRaw), num(p.pAdj), p.signif].join(","),
    ),
  ];
  downloadBlob(
    "qpcr_pairwise_tests.csv",
    new Blob(["\uFEFF" + lines.join("\n")], { type: "text/csv;charset=utf-8" }),
  );
}

/** 柱状图 PNG：SVG 序列化 → Image → canvas（2x 分辨率）→ toBlob 下载。 */
export async function exportChartPng(svg: SVGSVGElement, width: number, height: number): Promise<void> {
  const clone = svg.cloneNode(true) as SVGSVGElement;
  clone.setAttribute("xmlns", "http://www.w3.org/2000/svg");
  clone.setAttribute("width", String(width));
  clone.setAttribute("height", String(height));
  // 白底（PNG 需要不透明背景）
  const bg = document.createElementNS("http://www.w3.org/2000/svg", "rect");
  bg.setAttribute("x", "0");
  bg.setAttribute("y", "0");
  bg.setAttribute("width", String(width));
  bg.setAttribute("height", String(height));
  bg.setAttribute("fill", "#ffffff");
  clone.insertBefore(bg, clone.firstChild);

  const svgText = new XMLSerializer().serializeToString(clone);
  const svgUrl = URL.createObjectURL(new Blob([svgText], { type: "image/svg+xml;charset=utf-8" }));
  try {
    const img = await new Promise<HTMLImageElement>((resolve, reject) => {
      const el = new Image();
      el.onload = () => resolve(el);
      el.onerror = () => reject(new Error("SVG 渲染失败"));
      el.src = svgUrl;
    });
    const scale = 2;
    const canvas = document.createElement("canvas");
    canvas.width = width * scale;
    canvas.height = height * scale;
    const ctx = canvas.getContext("2d");
    if (!ctx) throw new Error("无法创建 canvas 上下文");
    ctx.fillStyle = "#ffffff";
    ctx.fillRect(0, 0, canvas.width, canvas.height);
    ctx.drawImage(img, 0, 0, canvas.width, canvas.height);
    const blob = await new Promise<Blob>((resolve, reject) => {
      canvas.toBlob((b) => (b ? resolve(b) : reject(new Error("PNG 编码失败"))), "image/png");
    });
    downloadBlob("qpcr_barplot.png", blob);
  } finally {
    URL.revokeObjectURL(svgUrl);
  }
}
