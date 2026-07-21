/**
 * Home.tsx — qPCR ΔΔCt 计算器主页面（design/kimi-ui-brief.md §10 页面结构）
 */
import { useMemo, useRef, useState } from "react";
import { Header } from "@/components/qpcr/Header.tsx";
import { InstrumentFlow } from "@/components/qpcr/InstrumentFlow.tsx";
import { BarChart, CHART_H, CHART_W } from "@/components/qpcr/BarChart.tsx";
import {
  Alert,
  Card,
  CardTitle,
  DataTable,
  FieldLabel,
  KimiButton,
  KimiSelect,
  Segmented,
  type TableCol,
} from "@/components/qpcr/ui.tsx";
import { readPastedTable, readTableFile } from "@/lib/fileio.ts";
import {
  ADJUST_LABELS,
  defaultColumns,
  fmt4,
  fmtP,
  runAnalysis,
  TEST_LABELS,
  type AnalysisResult,
  type PAdjust,
  type PreparedTable,
  type TestKind,
} from "@/lib/analysis.ts";
import { exportChartPng, exportCsv, exportExcel } from "@/lib/exporters.ts";

type InputMode = "paste" | "upload" | "instrument";

const MODE_OPTIONS: Array<{ value: InputMode; label: string }> = [
  { value: "paste", label: "复制粘贴" },
  { value: "upload", label: "上传整理表" },
  { value: "instrument", label: "仪器原始文件" },
];

const PASTE_PLACEHOLDER = "group\tGAPDH\tGENE\nNC\t12.1\t23.2\nNC\t12.3\t22.8\nTreatmentA\t12.4\t21.9";

export default function Home() {
  const [mode, setMode] = useState<InputMode>("paste");

  // ---- 三种输入来源各自的解析结果 ----
  const [pasteText, setPasteText] = useState("");
  const pasteParsed = useMemo((): { table?: PreparedTable; error?: string } => {
    if (!pasteText.trim()) return {};
    try {
      return { table: readPastedTable(pasteText) };
    } catch (e) {
      return { error: e instanceof Error ? e.message : String(e) };
    }
  }, [pasteText]);

  const [uploadTable, setUploadTable] = useState<PreparedTable | null>(null);
  const [uploadName, setUploadName] = useState("");
  const [uploadError, setUploadError] = useState("");

  const [instTable, setInstTable] = useState<PreparedTable | null>(null);

  const prepared: PreparedTable | null =
    mode === "paste" ? pasteParsed.table ?? null : mode === "upload" ? uploadTable : instTable;

  // ---- 分析设置（targetCol/hkCol 仅为用户覆盖值；生效值在渲染期派生，避免 setState-in-effect）----
  const [targetCol, setTargetCol] = useState("");
  const [hkCol, setHkCol] = useState("");
  const [ctrlGroup, setCtrlGroup] = useState("NC");
  const [testKind, setTestKind] = useState<TestKind>("welch");
  const [pAdjust, setPAdjust] = useState<PAdjust>("bh");
  const [yPadding, setYPadding] = useState(0.08);

  // ---- 计算结果 ----
  const [result, setResult] = useState<AnalysisResult | null>(null);
  const [calcError, setCalcError] = useState("");
  const svgRef = useRef<SVGSVGElement>(null);

  const numCandidates = useMemo(
    () => (prepared ? prepared.columns.filter((c) => c !== "group" && c !== "sample") : []),
    [prepared],
  );

  // 渲染期派生生效列：用户选择失效（列不存在或与管家列冲突）时回退到自动默认
  const defaults = useMemo(() => defaultColumns(numCandidates), [numCandidates]);
  const effHk = hkCol && numCandidates.includes(hkCol) ? hkCol : (defaults.hk ?? "");
  const effTarget =
    targetCol && numCandidates.includes(targetCol) && targetCol !== effHk
      ? targetCol
      : (defaults.target ?? "");

  const clearResult = () => {
    setResult(null);
    setCalcError("");
  };

  const onUploadFile = (f: File | undefined) => {
    if (!f) return;
    setUploadError("");
    clearResult();
    f.arrayBuffer()
      .then((buf) => {
        setUploadTable(readTableFile(f.name, buf));
        setUploadName(f.name);
      })
      .catch((e) => {
        setUploadTable(null);
        setUploadError(e instanceof Error ? e.message : String(e));
      });
  };

  const run = () => {
    setCalcError("");
    setResult(null);
    if (!prepared) return;
    try {
      setResult(
        runAnalysis(prepared, { targetCol: effTarget, hkCol: effHk, ctrlGroup: ctrlGroup.trim(), testKind, pAdjust, yPadding }),
      );
    } catch (e) {
      setCalcError(e instanceof Error ? e.message : String(e));
    }
  };

  const clearAll = () => {
    setPasteText("");
    setUploadTable(null);
    setUploadName("");
    setUploadError("");
    setInstTable(null);
    setResult(null);
    setCalcError("");
  };

  const hasGroup = prepared ? prepared.columns.includes("group") : false;

  // ---- 结果表列 ----
  const resultCols: TableCol[] = result
    ? result.columns.map((c) => ({
        key: c,
        title: c,
        numeric: ["dct", "ddct", "mrna"].includes(c),
        render: (row) => fmt4(row[c]),
      }))
    : [];

  const previewCols: TableCol[] = prepared
    ? prepared.columns.map((c) => ({ key: c, title: c, render: (row) => fmt4(row[c]) }))
    : [];

  return (
    <div className="min-h-screen bg-white">
      <Header />
      <main className="mx-auto max-w-[880px] px-4 pb-16 pt-24">
        {/* 标题区 */}
        <div className="mb-8">
          <h1 className="text-xl font-semibold leading-[30px] text-label-primary">
            qPCR ΔΔCt 计算器
          </h1>
          <p className="mt-1 text-sm leading-5 text-label-tertiary">
            从仪器原始文件或整理好的表格出发，完成 ΔΔCt 相对表达量计算、两两统计检验与作图，全部在浏览器本地运行。
          </p>
        </div>

        <div className="space-y-8">
          {/* 数据输入 */}
          <Card>
            <CardTitle hint="三种输入方式任选其一；仪器原始文件支持多通道多文件">
              数据输入
            </CardTitle>
            <div className="mb-5">
              <Segmented
                options={MODE_OPTIONS}
                value={mode}
                onChange={(v) => {
                  setMode(v);
                  clearResult();
                }}
              />
            </div>

            {mode === "paste" ? (
              <div>
                <textarea
                  value={pasteText}
                  onChange={(e) => {
                    setPasteText(e.target.value);
                    clearResult();
                  }}
                  placeholder={PASTE_PLACEHOLDER}
                  rows={7}
                  className="kimi-field w-full font-mono text-[13px] leading-5"
                />
                <p className="mt-2 text-xs leading-[18px] text-label-tertiary">
                  首行为列名，需包含 group 列 + 至少两列 Ct 数值（目标与管家），支持制表符 / 逗号分隔。
                </p>
                {pasteParsed.error ? (
                  <div className="mt-2">
                    <Alert kind="error">{pasteParsed.error}</Alert>
                  </div>
                ) : null}
              </div>
            ) : null}

            {mode === "upload" ? (
              <div>
                <label
                  className="flex cursor-pointer flex-col items-center justify-center rounded-lg border border-dashed border-sep-s1 bg-ground-pc px-4 py-6 transition-colors duration-150 hover:bg-fill-f1"
                >
                  <span className="text-sm leading-5 text-label-primary">
                    {uploadName ? `已选择：${uploadName}` : "点击选择整理好的 CSV / Excel 文件"}
                  </span>
                  <span className="mt-1 text-xs leading-[18px] text-label-tertiary">
                    结构同「复制粘贴」：group 列 + 各基因 Ct 列
                  </span>
                  <input
                    type="file"
                    accept=".csv,.tsv,.txt,.xls,.xlsx"
                    className="hidden"
                    onChange={(e) => {
                      onUploadFile(e.target.files?.[0]);
                      e.target.value = "";
                    }}
                  />
                </label>
                {uploadError ? (
                  <div className="mt-2">
                    <Alert kind="error">{uploadError}</Alert>
                  </div>
                ) : null}
              </div>
            ) : null}

            {/* 仪器模式保持挂载以保留子流程状态 */}
            <div className={mode === "instrument" ? "" : "hidden"}>
              <InstrumentFlow
                onPrepared={(table) => {
                  setInstTable(table);
                  clearResult();
                }}
              />
            </div>
          </Card>

          {/* 无数据提示 */}
          {!prepared && mode !== "instrument" ? (
            <Alert kind="info">
              请在上方粘贴数据或上传文件。至少包含：group 列 + 2 列 Ct 数值（目标与管家）。
            </Alert>
          ) : null}
          {!prepared && mode === "instrument" ? (
            <Alert kind="info">请完成上方仪器文件解析与分组步骤，点击「确认分组，生成分析表」继续。</Alert>
          ) : null}

          {/* 数据预览 */}
          {prepared ? (
            <Card className="result-enter">
              <div className="mb-4 flex items-center justify-between">
                <h2 className="text-base font-medium leading-6 text-label-primary">数据预览</h2>
                <KimiButton variant="dangerText" onClick={clearAll}>
                  清空
                </KimiButton>
              </div>
              <DataTable
                cols={previewCols}
                rows={prepared.rows.slice(0, 20) as Array<Record<string, unknown>>}
                maxHeight={360}
              />
              <p className="mt-2 text-xs leading-[18px] text-label-tertiary">
                共 {prepared.rows.length} 行 × {prepared.columns.length} 列（显示前 20 行）
              </p>
            </Card>
          ) : null}

          {/* 分析设置 */}
          {prepared && !hasGroup ? (
            <Alert kind="error">缺少列：group（区分分组）。请补充后重试。</Alert>
          ) : null}
          {prepared && hasGroup && numCandidates.length < 2 ? (
            <Alert kind="error">检测到的数值列不足（至少需要 2 列 Ct）。</Alert>
          ) : null}
          {prepared && hasGroup && numCandidates.length >= 2 ? (
            <Card className="result-enter">
              <CardTitle hint="目标基因与管家基因列、对照组、检验方式与作图参数">
                分析设置
              </CardTitle>
              <div className="grid grid-cols-1 gap-4 sm:grid-cols-2">
                <div>
                  <FieldLabel>目标基因 Ct 列</FieldLabel>
                  <KimiSelect value={effTarget} onChange={(e) => setTargetCol(e.target.value)}>
                    {numCandidates
                      .filter((c) => c !== effHk)
                      .map((c) => (
                        <option key={c} value={c}>{c}</option>
                      ))}
                  </KimiSelect>
                </div>
                <div>
                  <FieldLabel>管家基因 Ct 列</FieldLabel>
                  <KimiSelect value={effHk} onChange={(e) => setHkCol(e.target.value)}>
                    {numCandidates.map((c) => (
                      <option key={c} value={c}>{c}</option>
                    ))}
                  </KimiSelect>
                </div>
                <div>
                  <FieldLabel>对照组名称（与 group 列一致）</FieldLabel>
                  <input
                    value={ctrlGroup}
                    onChange={(e) => setCtrlGroup(e.target.value)}
                    className="kimi-field"
                    placeholder="NC"
                  />
                </div>
                <div>
                  <FieldLabel>两两比较检验方式</FieldLabel>
                  <KimiSelect value={testKind} onChange={(e) => setTestKind(e.target.value as TestKind)}>
                    {(Object.keys(TEST_LABELS) as TestKind[]).map((k) => (
                      <option key={k} value={k}>{TEST_LABELS[k]}</option>
                    ))}
                  </KimiSelect>
                </div>
                <div>
                  <FieldLabel>多重比较校正</FieldLabel>
                  <KimiSelect value={pAdjust} onChange={(e) => setPAdjust(e.target.value as PAdjust)}>
                    {(Object.keys(ADJUST_LABELS) as PAdjust[]).map((k) => (
                      <option key={k} value={k}>{ADJUST_LABELS[k]}</option>
                    ))}
                  </KimiSelect>
                </div>
                <div>
                  <FieldLabel>显著性横线与柱顶的间距（{yPadding.toFixed(2)}）</FieldLabel>
                  <input
                    type="range"
                    min={0.02}
                    max={0.3}
                    step={0.01}
                    value={yPadding}
                    onChange={(e) => setYPadding(Number(e.target.value))}
                    className="kimi-range mt-2"
                    aria-label="显著性横线与柱顶的间距"
                  />
                </div>
              </div>
            </Card>
          ) : null}

          {/* 主按钮 */}
          {prepared && hasGroup && numCandidates.length >= 2 ? (
            <div>
              <KimiButton variant="primary" size="hero" block onClick={run}>
                开始计算
              </KimiButton>
              {calcError ? (
                <div className="mt-3">
                  <Alert kind="error">{calcError}</Alert>
                </div>
              ) : null}
            </div>
          ) : null}

          {/* 结果区 */}
          {result ? (
            <>
              <div className="result-enter">
                <Alert kind="success">
                  计算完成：ΔCt = 目标 − 管家；ΔΔCt = ΔCt − 对照组 ΔCt 均值（nc_mean ={" "}
                  {fmt4(result.ncMean)}）；相对表达 = 2^(−ΔΔCt)
                </Alert>
              </div>

              <Card className="result-enter" >
                <CardTitle hint="展示保留 4 位小数；导出文件保留全精度">计算结果</CardTitle>
                <DataTable
                  cols={resultCols}
                  rows={result.rows as Array<Record<string, unknown>>}
                  maxHeight={420}
                />
              </Card>

              <Card className="result-enter">
                <CardTitle hint="按组名排序；SD 为样本标准差（ddof=1）">分组统计（mrna）</CardTitle>
                <DataTable
                  cols={[
                    { key: "group", title: "group" },
                    { key: "n", title: "N", numeric: true },
                    { key: "mean", title: "mrna_mean", numeric: true, render: (r) => fmt4(r.mean) },
                    { key: "sd", title: "mrna_sd", numeric: true, render: (r) => fmt4(r.sd) },
                  ]}
                  rows={result.gstats as unknown as Array<Record<string, unknown>>}
                />
              </Card>

              <Card className="result-enter">
                <CardTitle hint="横线标注显著性星号与校正后 p 值；* p<0.05，** p<0.01，*** p<0.001，**** p<0.0001">
                  作图
                </CardTitle>
                <div className="flex justify-center">
                  <div className="w-full max-w-[560px]">
                    <BarChart
                      ref={svgRef}
                      levels={result.levels}
                      means={result.means}
                      sds={result.sds}
                      pairwise={result.pairwise}
                      yPadding={result.settings.yPadding}
                    />
                  </div>
                </div>
                <p className="mt-3 text-center text-xs leading-[18px] text-label-tertiary">
                  柱高为各组相对表达量均值，误差线为 ±SD；{TEST_LABELS[result.settings.testKind]}，
                  {ADJUST_LABELS[result.settings.pAdjust]}。
                </p>
              </Card>

              <Card className="result-enter">
                <CardTitle hint="p_raw 为原始 p 值，p_adj 为多重校正后 p 值">
                  两两比较
                </CardTitle>
                <DataTable
                  cols={[
                    { key: "g1", title: "Group1" },
                    { key: "g2", title: "Group2" },
                    { key: "n1", title: "N1", numeric: true },
                    { key: "n2", title: "N2", numeric: true },
                    { key: "test", title: "Test" },
                    { key: "pRaw", title: "p_raw", numeric: true, render: (r) => fmtP(r.pRaw as number) },
                    { key: "pAdj", title: "p_adj", numeric: true, render: (r) => fmtP(r.pAdj as number) },
                    { key: "signif", title: "Signif" },
                  ]}
                  rows={result.pairwise as unknown as Array<Record<string, unknown>>}
                />
              </Card>

              <Card className="result-enter">
                <CardTitle hint="Excel 含 qpcr_results / settings / group_stats / pairwise_tests 四个 sheet">
                  导出
                </CardTitle>
                <div className="grid grid-cols-1 gap-3 sm:grid-cols-3">
                  <KimiButton variant="secondary" block onClick={() => exportExcel(result)}>
                    下载 Excel（4 个 sheet）
                  </KimiButton>
                  <KimiButton variant="secondary" block onClick={() => exportCsv(result)}>
                    下载两两比较 CSV
                  </KimiButton>
                  <KimiButton
                    variant="secondary"
                    block
                    onClick={() => {
                      if (svgRef.current) {
                        void exportChartPng(svgRef.current, CHART_W, CHART_H);
                      }
                    }}
                  >
                    下载柱状图 PNG（2x）
                  </KimiButton>
                </div>
              </Card>
            </>
          ) : null}
        </div>
      </main>

      <footer className="border-t border-sep-s1 bg-bg-secondary">
        <div className="mx-auto max-w-[880px] px-4 py-5 text-xs leading-[18px] text-label-tertiary">
          ΔCt = 目标基因 Ct − 管家基因 Ct；ΔΔCt = ΔCt − 对照组 ΔCt 均值；相对表达量 =
          2^(−ΔΔCt)。基于 Livak ΔΔCt 方法，全部计算在浏览器本地完成。
        </div>
      </footer>
    </div>
  );
}
