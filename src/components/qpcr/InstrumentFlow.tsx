/**
 * InstrumentFlow.tsx — 仪器原始文件模式子流程
 * 上传 → 逐文件解析报告 → （失败文件）手动解析 → 靶基因命名 →
 * 长表预览 → 分组映射（autoGroupSamples 三规则 + 可编辑）→
 * 管家基因候选提示 → 生成分析宽表
 */
import { useMemo, useRef, useState } from "react";
import {
  autoGroupSamples,
  dedupHeader,
  parseManualGrid,
  pivotToWide,
  targetsInAllSamples,
  InstrumentParseError,
  type FileReport,
  type Grid,
  type LongRow,
} from "@/lib/parsers.ts";
import { fileToGrid, parseInstrumentFile } from "@/lib/fileio.ts";
import { downloadInstrumentExample } from "@/lib/templates.ts";
import { fmt4, type PreparedTable } from "@/lib/analysis.ts";
import { Alert, CardTitle, DataTable, FieldLabel, KimiButton, KimiSelect } from "./ui.tsx";

interface InstFile {
  name: string;
  buf: ArrayBuffer;
}

type GroupRule = "identity" | "strip_suffix" | "prefix_sep";

const GROUP_RULES: Array<{ value: GroupRule; label: string }> = [
  { value: "identity", label: "每个样本独立成组" },
  { value: "strip_suffix", label: "前缀分组：去掉尾部重复编号（如 _1、-2、空格3）" },
  { value: "prefix_sep", label: "前缀分组：取第一个分隔符（_ - 空格）之前的部分" },
];

const NONE_OPT = "（无）";

// --------------------------------------------------------------------------
// 上传区（虚线 + 拖拽高亮）
// --------------------------------------------------------------------------
function UploadZone({ onFiles }: { onFiles: (files: File[]) => void }) {
  const [dragOver, setDragOver] = useState(false);
  const inputRef = useRef<HTMLInputElement>(null);
  return (
    <div
      role="button"
      tabIndex={0}
      onClick={() => inputRef.current?.click()}
      onKeyDown={(e) => {
        if (e.key === "Enter" || e.key === " ") inputRef.current?.click();
      }}
      onDragOver={(e) => {
        e.preventDefault();
        setDragOver(true);
      }}
      onDragLeave={() => setDragOver(false)}
      onDrop={(e) => {
        e.preventDefault();
        setDragOver(false);
        onFiles(Array.from(e.dataTransfer.files));
      }}
      className={`flex cursor-pointer flex-col items-center justify-center rounded-lg border border-dashed px-4 py-8 transition-colors duration-150 ${
        dragOver ? "border-kimi-blue bg-light-blue-bg" : "border-sep-s1 bg-ground-pc hover:bg-fill-f1"
      }`}
    >
      <p className="text-sm leading-5 text-label-primary">点击选择或拖拽仪器导出文件到此处</p>
      <p className="mt-1 text-xs leading-[18px] text-label-tertiary">
        支持多选：LightCycler 480 / Rotor-Gene 等分通道导出的仪器请一次上传全部通道文件
      </p>
      <p className="mt-1 text-xs leading-[18px] text-label-tertiary">
        .eds / .sds / .pcrd 为二进制工程文件，请先在仪器软件中 Export 为 txt / csv / xlsx
      </p>
      <input
        ref={inputRef}
        type="file"
        multiple
        accept=".txt,.csv,.tsv,.xls,.xlsx,.eds,.sds,.pcrd"
        className="hidden"
        onChange={(e) => {
          onFiles(Array.from(e.target.files ?? []));
          e.target.value = "";
        }}
      />
    </div>
  );
}

// --------------------------------------------------------------------------
// 手动解析面板（自动识别失败的文件：选表头行 + well/sample/target/ct 列）
// --------------------------------------------------------------------------
function ManualParsePanel({
  file,
  onParsed,
}: {
  file: InstFile;
  onParsed: (rows: LongRow[], report: FileReport) => void;
}) {
  const [headerRow, setHeaderRow] = useState(1); // 1 基
  // 列选择：以「列集合」为 key 的袋式状态；表头行变化导致列集合变化时渲染期自动回退默认，
  // 无需 useEffect 重置（避免 setState-in-effect 循环风险）
  const [selBag, setSelBag] = useState<{
    key: string;
    well: string;
    sample: string;
    target: string;
    ct: string | null;
  }>({ key: "", well: NONE_OPT, sample: NONE_OPT, target: NONE_OPT, ct: null });
  const [tname, setTname] = useState("");
  const [error, setError] = useState("");

  const gridOrError = useMemo((): { grid?: Grid; error?: string } => {
    try {
      return { grid: fileToGrid(file.name, file.buf) };
    } catch (e) {
      return { error: e instanceof Error ? e.message : String(e) };
    }
  }, [file]);

  const grid = gridOrError.grid;
  const cols = useMemo(
    () => (grid && headerRow >= 1 && headerRow <= grid.length ? dedupHeader(grid[headerRow - 1]) : []),
    [grid, headerRow],
  );
  const colsKey = cols.join("");
  const sel = selBag.key === colsKey ? selBag : null;
  const wellCol = sel && cols.includes(sel.well) ? sel.well : NONE_OPT;
  const sampleCol = sel && cols.includes(sel.sample) ? sel.sample : NONE_OPT;
  const targetCol = sel && cols.includes(sel.target) ? sel.target : NONE_OPT;
  const ctCol = sel && sel.ct !== null && cols.includes(sel.ct) ? sel.ct : (cols[cols.length - 1] ?? "");
  const updateSel = (patch: Partial<{ well: string; sample: string; target: string; ct: string }>) =>
    setSelBag({ key: colsKey, well: wellCol, sample: sampleCol, target: targetCol, ct: ctCol, ...patch });

  if (gridOrError.error) {
    return <Alert kind="error">文件无法读成表格：{gridOrError.error}</Alert>;
  }
  if (!grid) return null;

  const lo = Math.max(1, headerRow - 2);
  const hi = Math.min(grid.length, headerRow + 5);
  const previewRows = grid.slice(lo - 1, hi).map((row, i) => {
    const rec: Record<string, unknown> = { __no: lo + i };
    row.forEach((c, j) => {
      rec[`c${j}`] = c;
    });
    return rec;
  });
  const maxCols = Math.max(1, ...grid.slice(lo - 1, hi).map((r) => r.length));
  const previewCols = [
    { key: "__no", title: "行号" },
    ...Array.from({ length: maxCols }, (_, j) => ({ key: `c${j}`, title: `列${j + 1}` })),
  ];

  const doParse = () => {
    setError("");
    try {
      const { rows, report } = parseManualGrid(
        grid,
        file.name,
        headerRow - 1,
        {
          well: wellCol === NONE_OPT ? undefined : wellCol,
          sample: sampleCol === NONE_OPT ? undefined : sampleCol,
          target: targetCol === NONE_OPT ? undefined : targetCol,
          ct: ctCol || undefined,
        },
        tname.trim() || undefined,
      );
      onParsed(rows, report);
    } catch (e) {
      setError(e instanceof InstrumentParseError ? e.message : String(e));
    }
  };

  return (
    <div className="space-y-3">
      <div className="w-48">
        <FieldLabel>表头所在行（从 1 开始，共 {grid.length} 行）</FieldLabel>
        <input
          type="number"
          min={1}
          max={grid.length}
          value={headerRow}
          onChange={(e) => setHeaderRow(Math.max(1, Math.min(grid.length, Number(e.target.value) || 1)))}
          className="kimi-field"
        />
      </div>
      <div>
        <p className="mb-1.5 text-xs leading-[18px] text-label-tertiary">
          第 {lo}–{hi} 行预览：
        </p>
        <DataTable cols={previewCols} rows={previewRows} maxHeight={200} />
      </div>
      <div className="grid grid-cols-2 gap-3 sm:grid-cols-4">
        <div>
          <FieldLabel optional>孔位列</FieldLabel>
          <KimiSelect value={wellCol} onChange={(e) => updateSel({ well: e.target.value })}>
            {[NONE_OPT, ...cols].map((c) => (
              <option key={c} value={c}>{c}</option>
            ))}
          </KimiSelect>
        </div>
        <div>
          <FieldLabel optional>样本列</FieldLabel>
          <KimiSelect value={sampleCol} onChange={(e) => updateSel({ sample: e.target.value })}>
            {[NONE_OPT, ...cols].map((c) => (
              <option key={c} value={c}>{c}</option>
            ))}
          </KimiSelect>
        </div>
        <div>
          <FieldLabel optional>靶基因列</FieldLabel>
          <KimiSelect value={targetCol} onChange={(e) => updateSel({ target: e.target.value })}>
            {[NONE_OPT, ...cols].map((c) => (
              <option key={c} value={c}>{c}</option>
            ))}
          </KimiSelect>
        </div>
        <div>
          <FieldLabel>Ct/Cq 列</FieldLabel>
          <KimiSelect value={ctCol} onChange={(e) => updateSel({ ct: e.target.value })}>
            {cols.map((c) => (
              <option key={c} value={c}>{c}</option>
            ))}
          </KimiSelect>
        </div>
      </div>
      {targetCol === NONE_OPT ? (
        <div>
          <FieldLabel>靶基因名称（无靶基因列时填写）</FieldLabel>
          <input
            value={tname}
            onChange={(e) => setTname(e.target.value)}
            placeholder="例如 GAPDH"
            className="kimi-field"
          />
        </div>
      ) : null}
      {error ? <Alert kind="error">{error}</Alert> : null}
      <KimiButton variant="primary" block onClick={doParse}>
        解析并加入数据
      </KimiButton>
    </div>
  );
}

// --------------------------------------------------------------------------
// 单选（Kimi 风格 16px 圆点）
// --------------------------------------------------------------------------
function RadioRow({
  checked,
  onChange,
  label,
}: {
  checked: boolean;
  onChange: () => void;
  label: string;
}) {
  return (
    <label className="flex cursor-pointer items-center gap-2 py-1 text-sm leading-5 text-label-primary">
      <input type="radio" checked={checked} onChange={onChange} className="sr-only" />
      <span
        aria-hidden
        className={`flex h-4 w-4 items-center justify-center rounded-full border transition-colors duration-150 ${
          checked ? "border-kimi-blue" : "border-sep-s1"
        }`}
      >
        {checked ? <span className="h-2 w-2 rounded-full bg-kimi-blue" /> : null}
      </span>
      {label}
    </label>
  );
}

// --------------------------------------------------------------------------
// 主流程
// --------------------------------------------------------------------------
export function InstrumentFlow({
  onPrepared,
}: {
  onPrepared: (table: PreparedTable, sourceDesc: string) => void;
}) {
  const [files, setFiles] = useState<InstFile[]>([]);
  const [rows, setRows] = useState<LongRow[]>([]);
  const [reports, setReports] = useState<FileReport[]>([]);
  const [renames, setRenames] = useState<Record<string, string>>({});
  const [groupRule, setGroupRule] = useState<GroupRule>("identity");
  // 分组手动编辑：以「样本集+规则」为 key 的袋式状态，key 变化时渲染期自动失效，
  // 无需 useEffect 重置（避免 setState-in-effect 循环风险）
  const [editBag, setEditBag] = useState<{ key: string; edits: Record<string, string> }>({
    key: "",
    edits: {},
  });
  const [previewOpen, setPreviewOpen] = useState(false);
  const [manualFile, setManualFile] = useState("");
  const [generatedMsg, setGeneratedMsg] = useState("");
  const [parseError, setParseError] = useState("");

  // 选择文件后立即解析全部
  const addFiles = (fl: File[]) => {
    if (fl.length === 0) return;
    setParseError("");
    setGeneratedMsg("");
    Promise.all(
      fl.map(async (f) => ({ name: f.name, buf: await f.arrayBuffer() })),
    )
      .then((instFiles) => {
        setFiles(instFiles);
        const allRows: LongRow[] = [];
        const allReports: FileReport[] = [];
        const initRenames: Record<string, string> = {};
        for (const f of instFiles) {
          try {
            const { rows: r, report } = parseInstrumentFile(f.name, f.buf);
            allRows.push(...r);
            allReports.push(report);
            if (report.needs_target_name) {
              initRenames[f.name] = report.default_target || f.name;
            }
          } catch (e) {
            allReports.push({
              file: f.name,
              instrument: "未知",
              status: "error",
              rows: 0,
              targets: [],
              needs_target_name: false,
              default_target: "",
              message: e instanceof Error ? e.message : String(e),
            });
          }
        }
        setRows(allRows);
        setReports(allReports);
        setRenames(initRenames);
        setEditBag({ key: "", edits: {} });
        setManualFile(allReports.find((r) => r.status === "error")?.file ?? "");
      })
      .catch((e) => setParseError(String(e)));
  };

  // 应用靶基因命名后的有效长表
  const effectiveRows = useMemo(() => {
    const needs = new Set(reports.filter((r) => r.needs_target_name).map((r) => r.file));
    if (needs.size === 0) return rows;
    return rows.map((r) =>
      needs.has(r.source_file) && renames[r.source_file]
        ? { ...r, target: renames[r.source_file] }
        : r,
    );
  }, [rows, reports, renames]);

  const samples = useMemo(
    () => [...new Set(effectiveRows.map((r) => String(r.sample)))].sort(),
    [effectiveRows],
  );
  const samplesKey = samples.join("");
  const baseMap = useMemo(() => autoGroupSamples(samples, groupRule), [samplesKey, groupRule]); // eslint-disable-line react-hooks/exhaustive-deps
  const groupKey = `${samplesKey}|${groupRule}`;
  const groupEdits = editBag.key === groupKey ? editBag.edits : {};
  const groupMap: Record<string, string> = { ...baseMap, ...groupEdits };
  const setGroupEdit = (sample: string, value: string) =>
    setEditBag({ key: groupKey, edits: { ...groupEdits, [sample]: value } });

  const failedReports = reports.filter((r) => r.status === "error");
  const failedFiles = files.filter((f) => failedReports.some((r) => r.file === f.name));
  const needsNaming = reports.filter((r) => r.needs_target_name);

  // 管家基因候选：出现于所有样本的靶标
  const sortedTargets = useMemo(() => targetsInAllSamples(effectiveRows), [effectiveRows]);
  const presentAll = useMemo(() => {
    const allS = new Set(samples);
    return sortedTargets.filter(
      (t) =>
        new Set(effectiveRows.filter((r) => r.target === t).map((r) => String(r.sample))).size ===
        allS.size,
    );
  }, [sortedTargets, effectiveRows, samples]);

  const generate = () => {
    const wide = pivotToWide(effectiveRows.map((r) => ({ ...r, sample: String(r.sample) })));
    const targets = wide.columns.slice(1);
    const outRows = wide.rows.map((rec) => {
      const s = String(rec.sample);
      const g = (groupMap[s] ?? s).trim() || s;
      return { sample: s, group: g, ...Object.fromEntries(targets.map((t) => [t, rec[t]])) };
    });
    const table: PreparedTable = { columns: ["sample", "group", ...targets], rows: outRows };
    setGeneratedMsg(`分析表已生成：${outRows.length} 个样本 × ${targets.length} 个基因列。`);
    onPrepared(table, "仪器原始文件");
  };

  const longPreviewCols = [
    { key: "well", title: "孔位" },
    { key: "sample", title: "样本" },
    { key: "target", title: "靶基因" },
    { key: "ct", title: "Ct", numeric: true, render: (r: Record<string, unknown>) => fmt4(r.ct) },
    { key: "source_file", title: "来源文件" },
  ];

  const groupTableCols = [
    { key: "sample", title: "样本" },
    {
      key: "group",
      title: "分组（可编辑）",
      render: (r: Record<string, unknown>) => {
        const s = String(r.sample);
        return (
          <input
            value={groupMap[s] ?? s}
            onChange={(e) => setGroupEdit(s, e.target.value)}
            className="h-7 w-full min-w-[120px] rounded-md border border-sep-s1 bg-white px-2 text-sm outline-none transition-[border-color,box-shadow] duration-150 focus-visible:border-kimi-blue focus-visible:shadow-[0_0_0_1px_#1783ff]"
          />
        );
      },
    },
  ];

  return (
    <div className="space-y-5">
      <UploadZone onFiles={addFiles} />
      <div className="flex flex-wrap items-center gap-2">
        <span className="text-xs leading-[18px] text-label-tertiary">
          不确定该传什么文件？下载示例看看：
        </span>
        <KimiButton variant="outline" onClick={() => downloadInstrumentExample("quantstudio")}>
          QuantStudio 示例（.txt）
        </KimiButton>
        <KimiButton variant="outline" onClick={() => downloadInstrumentExample("biorad")}>
          Bio-Rad CFX 示例（.csv）
        </KimiButton>
        <span className="text-xs leading-[18px] text-label-tertiary">
          可直接上传到上方体验完整流程
        </span>
      </div>
      {parseError ? <Alert kind="error">{parseError}</Alert> : null}

      {reports.length > 0 ? (
        <div className="space-y-2">
          <CardTitle hint="每个文件的识别仪器、解析行数；失败文件附中文导出指引">解析报告</CardTitle>
          {reports.map((rep) =>
            rep.status === "ok" ? (
              <Alert key={rep.file} kind="success">
                <strong>{rep.file}</strong> — {rep.message}
                {rep.needs_target_name ? "（未识别到靶基因列，请在下方命名）" : ""}
              </Alert>
            ) : (
              <Alert key={rep.file} kind="error">
                <strong>{rep.file}</strong> — {rep.message}
              </Alert>
            ),
          )}
        </div>
      ) : null}

      {failedFiles.length > 0 ? (
        <div className="rounded-md border border-sep-s1 p-4">
          <CardTitle hint="自动识别失败的文件：自选表头行与 well / sample / target / ct 列">
            手动解析
          </CardTitle>
          <div className="mb-3 w-64">
            <FieldLabel>选择文件</FieldLabel>
            <KimiSelect value={manualFile} onChange={(e) => setManualFile(e.target.value)}>
              {failedFiles.map((f) => (
                <option key={f.name} value={f.name}>{f.name}</option>
              ))}
            </KimiSelect>
          </div>
          {manualFile && failedFiles.some((f) => f.name === manualFile) ? (
            <ManualParsePanel
              key={manualFile}
              file={failedFiles.find((f) => f.name === manualFile)!}
              onParsed={(newRows, rep) => {
                setRows((prev) => [...prev.filter((r) => r.source_file !== rep.file), ...newRows]);
                setReports((prev) => prev.map((r) => (r.file === rep.file ? rep : r)));
                // 解析成功后切到下一个仍失败的文件；没有则收起面板
                const nextReports = reports.map((r) => (r.file === rep.file ? rep : r));
                setManualFile(nextReports.find((r) => r.status === "error")?.file ?? "");
                if (rep.needs_target_name) {
                  setRenames((prev) => ({
                    ...prev,
                    [rep.file]: rep.default_target || prev[rep.file] || rep.file,
                  }));
                }
              }}
            />
          ) : null}
        </div>
      ) : null}

      {needsNaming.length > 0 ? (
        <div className="rounded-md border border-sep-s1 p-4">
          <CardTitle hint="以下文件按通道单独导出、文件内没有靶基因列，请为每个文件填写对应的靶基因名称">
            靶基因命名
          </CardTitle>
          <div className="grid grid-cols-1 gap-3 sm:grid-cols-2">
            {needsNaming.map((rep) => (
              <div key={rep.file}>
                <FieldLabel>「{rep.file}」的靶基因名称</FieldLabel>
                <input
                  value={renames[rep.file] ?? ""}
                  onChange={(e) =>
                    setRenames((prev) => ({ ...prev, [rep.file]: e.target.value }))
                  }
                  className="kimi-field"
                />
              </div>
            ))}
          </div>
        </div>
      ) : null}

      {effectiveRows.length > 0 ? (
        <div className="space-y-4">
          <div>
            <button
              type="button"
              onClick={() => setPreviewOpen((v) => !v)}
              className="flex items-center gap-1 text-sm leading-5 text-label-secondary transition-colors duration-150 hover:text-label-primary"
            >
              <span
                aria-hidden
                className={`inline-block transition-transform duration-150 ${previewOpen ? "rotate-90" : ""}`}
              >
                ▶
              </span>
              长表预览（前 30 行）
            </button>
            <p className="mt-1 text-xs leading-[18px] text-label-tertiary">
              共 {effectiveRows.length} 行 · {samples.length} 个样本 ·{" "}
              {new Set(effectiveRows.map((r) => r.target)).size} 个靶标（未检出 Ct 显示为空）
            </p>
            {previewOpen ? (
              <div className="mt-2">
                <DataTable
                  cols={longPreviewCols}
                  rows={effectiveRows.slice(0, 30) as unknown as Array<Record<string, unknown>>}
                  maxHeight={320}
                />
              </div>
            ) : null}
          </div>

          <div>
            <CardTitle hint="自动分组规则（生成后可逐行修改）">样本分组映射</CardTitle>
            <div className="mb-3 space-y-0.5">
              {GROUP_RULES.map((r) => (
                <RadioRow
                  key={r.value}
                  checked={groupRule === r.value}
                  onChange={() => setGroupRule(r.value)}
                  label={r.label}
                />
              ))}
            </div>
            <DataTable
              cols={groupTableCols}
              rows={samples.map((s) => ({ sample: s }))}
              maxHeight={280}
            />
          </div>

          {effectiveRows.length > 0 ? (
            presentAll.length > 0 ? (
              <Alert kind="info">
                出现于所有样本的靶标（适合作管家基因）：<strong>{presentAll.join("、")}</strong>
              </Alert>
            ) : (
              <Alert kind="warning">
                没有靶标覆盖全部样本，请检查是否漏传了某个通道的文件。
              </Alert>
            )
          ) : null}

          <KimiButton variant="primary" block onClick={generate}>
            确认分组，生成分析表
          </KimiButton>
          {generatedMsg ? <Alert kind="success">{generatedMsg}</Alert> : null}
        </div>
      ) : reports.length > 0 ? (
        <Alert kind="warning">尚未解析出有效数据行，请检查上方解析报告。</Alert>
      ) : null}
    </div>
  );
}
