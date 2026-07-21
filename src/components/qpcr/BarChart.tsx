/**
 * BarChart.tsx — 自绘 SVG 分组均值±SD 柱状图 + 两两显著性横线
 *
 * 横线分层避让算法移植自 qPCR-app/app.py 第 445–470 行 find_free_y：
 * 同一 x 区间内横线垂直间距不足 bump*0.9 时逐层上移。
 * 配色：同色相不透明度序列（design/kimi-ui-brief.md §8）。
 * 标题由外层卡片提供，SVG 内不画标题；顶部留白按最上层横线 + 两行标注自适应。
 */
import { forwardRef } from "react";
import { fmtG3, type PairRow } from "@/lib/analysis.ts";

export const CHART_W = 640;
export const CHART_H = 440;

const PALETTE = ["#1783ff", "#5ca8ff", "#8bc1ff", "#a2cdff", "#c5e0ff"];
const AXIS_TEXT = "rgba(0,0,0,0.6)";
const TICK_TEXT = "rgba(0,0,0,0.45)";
const GRID = "rgba(0,0,0,0.13)";
const BAR_EDGE = "rgba(0,0,0,0.75)";
const FONT = '-apple-system, "PingFang SC", "Helvetica Neue", "Microsoft YaHei", sans-serif';

/** 1/2/2.5/5 × 10^k 的 nice step。 */
function niceStep(raw: number): number {
  const mag = Math.pow(10, Math.floor(Math.log10(raw)));
  const norm = raw / mag;
  if (norm <= 1) return mag;
  if (norm <= 2) return 2 * mag;
  if (norm <= 2.5) return 2.5 * mag;
  if (norm <= 5) return 5 * mag;
  return 10 * mag;
}

interface BarChartProps {
  levels: string[];
  means: number[];
  sds: number[];
  pairwise: PairRow[];
  /** 显著性横线与柱顶的间距（相对比例） */
  yPadding: number;
}

export const BarChart = forwardRef<SVGSVGElement, BarChartProps>(function BarChart(
  { levels, means, sds, pairwise, yPadding },
  ref,
) {
  const ml = 72;
  const mr = 20;
  const mt = 24;
  const mb = 56;
  const plotW = CHART_W - ml - mr;
  const plotH = CHART_H - mt - mb;
  const n = levels.length;

  // ---- 数据范围（NaN 安全：全未检出组不参与纵轴上限）----
  // 柱顶+SD（仅有限值），用于 ymax_data
  const topsFinite = means.map((m, i) =>
    Number.isFinite(m) ? m + (Number.isFinite(sds[i]) ? sds[i] : 0) : NaN,
  );
  // 横线基点用：非有限 mean 按 0 计（对齐旧版 nan 不参与但不崩溃的行为）
  const tops0 = means.map(
    (m, i) => (Number.isFinite(m) ? m : 0) + (Number.isFinite(sds[i]) ? sds[i] : 0),
  );
  let ymaxData = Math.max(...topsFinite.filter(Number.isFinite));
  if (!Number.isFinite(ymaxData)) {
    const fm = means.filter(Number.isFinite);
    ymaxData = fm.length > 0 ? Math.max(...fm) : 1;
  }
  if (ymaxData <= 0) ymaxData = 1;
  const yMaxPre = ymaxData * (1 + Math.max(0.25, yPadding * pairwise.length));

  // ---- find_free_y 移植（x 以组索引为单位，y 以数值为单位；与坐标缩放无关）----
  const usedSpans: Array<[number, number, number]> = [];
  const findFreeY = (xl: number, xr: number, base: number): number => {
    let h = base;
    const bump = ymaxData * yPadding;
    while (
      usedSpans.some(([uxl, uxr, uy]) => xl <= uxr && xr >= uxl && Math.abs(h - uy) < bump * 0.9)
    ) {
      h += bump;
    }
    usedSpans.push([xl, xr, h]);
    return h;
  };

  const brackets = pairwise.map((p) => {
    const j1 = levels.indexOf(p.g1);
    const j2 = levels.indexOf(p.g2);
    const x1 = Math.min(j1, j2);
    const x2 = Math.max(j1, j2);
    const topY = Math.max(tops0[x1], tops0[x2]);
    const baseY = topY + ymaxData * yPadding * 0.4;
    const y = findFreeY(x1 - 0.2, x2 + 0.2, baseY);
    const dy = ymaxData * 0.01;
    const sub = Number.isNaN(p.pRaw) ? "" : `(p=${fmtG3(p.pAdj)})`;
    return { x1, x2, y, dy, label: p.signif, sub, key: `${p.g1}|${p.g2}` };
  });

  // ---- 顶部留白自适应：保证最上层横线 + 星号/p 值两行文字完整在画布内 ----
  const lineH = (16 / plotH) * yMaxPre; // 一行 12px 文字（含行距）折算成数据单位
  const requiredTop = brackets.reduce(
    (acc, b) => Math.max(acc, b.y + b.dy + ymaxData * 0.012 + (b.sub ? 2.3 : 1.2) * lineH),
    0,
  );
  const yMax = Math.max(yMaxPre, requiredTop);

  const step = niceStep(yMax / 4);
  const ticks: number[] = [];
  for (let v = 0; v <= yMax + step * 1e-9; v += step) ticks.push(v);
  const axisMax = ticks[ticks.length - 1];

  const X = (i: number): number => ml + (plotW / n) * (i + 0.5);
  const Y = (v: number): number => mt + plotH - (v / axisMax) * plotH;
  const barW = (plotW / n) * 0.6;

  const fmtTick = (v: number): string => String(Math.round(v * 1e6) / 1e6);

  return (
    <svg
      ref={ref}
      viewBox={`0 0 ${CHART_W} ${CHART_H}`}
      width={CHART_W}
      height={CHART_H}
      role="img"
      aria-label="qPCR 分组比较柱状图"
      style={{ maxWidth: "100%", height: "auto", fontFamily: FONT }}
      fontFamily={FONT}
    >
      {/* 横向虚线网格 + y 轴刻度 */}
      {ticks.map((v) => (
        <g key={v}>
          {v > 0 ? (
            <line
              x1={ml}
              x2={CHART_W - mr}
              y1={Y(v)}
              y2={Y(v)}
              stroke={GRID}
              strokeWidth={0.5}
              strokeDasharray="4 4"
            />
          ) : null}
          <text x={ml - 8} y={Y(v) + 4} textAnchor="end" fontSize={12} fill={TICK_TEXT}>
            {fmtTick(v)}
          </text>
        </g>
      ))}

      {/* 坐标轴（仅左/下） */}
      <line x1={ml} x2={ml} y1={mt} y2={mt + plotH} stroke={AXIS_TEXT} strokeWidth={0.75} />
      <line x1={ml} x2={CHART_W - mr} y1={mt + plotH} y2={mt + plotH} stroke={AXIS_TEXT} strokeWidth={0.75} />
      <text
        x={18}
        y={mt + plotH / 2}
        textAnchor="middle"
        fontSize={12}
        fill={AXIS_TEXT}
        transform={`rotate(-90 18 ${mt + plotH / 2})`}
      >
        Relative mRNA Expression
      </text>

      {/* 柱 + 误差线（均值非有限的组只画组名） */}
      {levels.map((lv, i) => {
        const m = means[i];
        const sd = Number.isFinite(sds[i]) ? sds[i] : NaN;
        const cx = X(i);
        const hasBar = Number.isFinite(m) && m > 0;
        const yTop = hasBar ? Y(m) : mt + plotH;
        const h = mt + plotH - yTop;
        return (
          <g key={lv}>
            {hasBar ? (
              <rect
                x={cx - barW / 2}
                y={yTop}
                width={barW}
                height={h}
                fill={PALETTE[i % PALETTE.length]}
                fillOpacity={0.92}
                stroke={BAR_EDGE}
                strokeWidth={0.6}
              />
            ) : null}
            {hasBar && Number.isFinite(sd) ? (
              <g stroke={BAR_EDGE} strokeWidth={1}>
                <line x1={cx} x2={cx} y1={Y(Math.max(0, m - sd))} y2={Y(m + sd)} />
                <line x1={cx - 5} x2={cx + 5} y1={Y(m + sd)} y2={Y(m + sd)} />
                <line x1={cx - 5} x2={cx + 5} y1={Y(Math.max(0, m - sd))} y2={Y(Math.max(0, m - sd))} />
              </g>
            ) : null}
            <text x={cx} y={mt + plotH + 20} textAnchor="middle" fontSize={12} fill={AXIS_TEXT}>
              {lv}
            </text>
          </g>
        );
      })}

      {/* 显著性横线 + 星号 / 校正后 p 值 */}
      {brackets.map((b) => {
        const px1 = X(b.x1);
        const px2 = X(b.x2);
        const yTop = Y(b.y + b.dy);
        const yBot = Y(b.y);
        const midX = (px1 + px2) / 2;
        const labelY = Y(b.y + ymaxData * 0.012) - 2;
        return (
          <g key={b.key}>
            <polyline
              points={`${px1},${yBot} ${px1},${yTop} ${px2},${yTop} ${px2},${yBot}`}
              fill="none"
              stroke={BAR_EDGE}
              strokeWidth={1}
            />
            <text x={midX} y={labelY} textAnchor="middle" fontSize={12} fill={AXIS_TEXT}>
              {b.sub ? (
                <>
                  <tspan x={midX} dy="-1.1em" fontWeight={500}>
                    {b.label}
                  </tspan>
                  <tspan x={midX} dy="1.25em">
                    {b.sub}
                  </tspan>
                </>
              ) : (
                b.label
              )}
            </text>
          </g>
        );
      })}
    </svg>
  );
});
