/**
 * parsers.ts — qPCR 仪器原始导出文件解析模块（浏览器端 TypeScript 移植版）
 *
 * 移植自 qPCR-app/parsers.py（Streamlit 版），行为逐条对齐：
 *   文本 → 二维网格(grid) → 表头行扫描 → 仪器格式识别 →
 *   列名同义词映射 → 过滤 NTC/Omit → 标准化长表。
 *
 * 与 Python 版的接口差异（架构要求）：
 *   - 核心函数接收「二维字符串网格 + 文件名」，不依赖浏览器 File / pandas；
 *     Excel 解码由应用层（SheetJS）完成后同样喂入网格函数。
 *   - ct 未检出用 null 表示（对应 Python 版的 NaN）。
 *   - parseTextFile(filename, text) 负责 txt/csv：转码由调用方完成
 *     （浏览器端 TextDecoder；可先用导出的 decodeText 做 utf-8→gb18030 嗅探），
 *     本函数内做分隔符嗅探、欧区小数逗号清洗、HTML 伪 xls 表格提取。
 *
 * 仅使用「可擦除」TypeScript 语法（无 enum / namespace / 参数属性），
 * 可被 Node 24 原生 type-stripping 直接执行。
 *
 * 支持仪器：
 *   - ABI/Thermo QuantStudio（QS 3/5/6/7/12K Flex、ViiA 7、7500、StepOne）v1.x / v2.x
 *   - Bio-Rad CFX（CFX96/384/Opus/Connect/Duet，CFX Maestro，含欧区分号+小数逗号）
 *   - Roche LightCycler 480（每通道一个文件）/ LightCycler 96（单文件含 Gene 列）
 *   - Agilent AriaMx（多 sheet Excel，应用层选 sheet 后喂网格）
 *   - Qiagen Rotor-Gene Q（每通道单独导出）
 *   - Eppendorf realplex（同义词兼容）
 *   - 国产仪器（天隆 Gentier、博日 LineGene、宏石 SLAN，中文表头）
 */

// --------------------------------------------------------------------------
// 类型
// --------------------------------------------------------------------------
export type Grid = string[][];

export type Field = "well" | "sample" | "target" | "ct";

/** 标准化长表行。ct 为 number，未检出为 null（Python 版为 NaN）。 */
export interface LongRow {
  well: string;
  sample: string;
  target: string;
  ct: number | null;
  source_file: string;
}

export interface FileReport {
  file: string;
  instrument: string;
  status: "ok" | "error";
  rows: number;
  targets: string[];
  needs_target_name: boolean;
  default_target: string;
  message: string;
}

export interface ParseResult {
  rows: LongRow[];
  report: FileReport;
}

export interface MultiParseResult {
  rows: LongRow[];
  reports: FileReport[];
}

/** _grid_to_raw_df 的产物：去重后的表头 + 数据行（全部为字符串）。 */
export interface RawTable {
  columns: string[];
  rows: string[][];
}

export interface ManualColumns {
  well?: string;
  sample?: string;
  target?: string;
  ct?: string;
}

// --------------------------------------------------------------------------
// 异常
// --------------------------------------------------------------------------
export class InstrumentParseError extends Error {
  constructor(message: string) {
    super(message);
    this.name = "InstrumentParseError";
  }
}

// --------------------------------------------------------------------------
// 列名同义词表（中英文，按优先级排序，匹配时做归一化后的精确匹配）
// 注意：与 parsers.py 保持逐字符一致（"cт" 第二个字符是西里尔字母 U+0442）。
// --------------------------------------------------------------------------
export const SYNONYMS: Record<Field, string[]> = {
  well: [
    "well", "well position", "wellposition", "well pos", "pos", "position",
    "well no", "well number", "no.", "no", "编号", "孔", "孔位", "孔号", "孔位置",
  ],
  sample: [
    "sample name", "samplename", "sample", "name", "sample id",
    "样本名称", "样品名称", "样本", "样品", "名称", "样本名", "样品名",
  ],
  target: [
    "target name", "targetname", "target", "gene", "detector name",
    "detector", "fluorophore", "fluor", "dye", "reporter", "channel",
    "靶基因", "目的基因", "基因", "检测基因", "荧光通道", "通道", "荧光", "染料",
  ],
  ct: [
    "ct", "cт", "cq", "cp", "crossing point", "crossingpoint",
    "ct值", "cq值", "cp值", "ct 值", "cq 值", "阈值循环数",
  ],
};

/** ct 兜底：只有均值列时的次选（优先级低） */
export const CT_FALLBACK: string[] = ["ct mean", "cq mean", "cp mean", "平均ct", "平均 ct", "平均cq"];

/** 未检出标记 → null */
export const UNDETERMINED: ReadonlySet<string> = new Set([
  "", "undetermined", "undeterm", "n/a", "na", "n.a.", "negative", "neg",
  "-", "--", "—", "none", "nan", "null", "no ct", "no cq", "no cp",
  "not detected", "no amplification", "no amp", "undet", "未检出", "无",
]);

/** 二进制工程文件（不可解析，需引导用户先在仪器软件里 Export） */
export const BINARY_HINTS: Record<string, string> = {
  eds: "这是 ABI QuantStudio 的实验工程文件（.eds），无法直接解析。请在 QuantStudio 软件中打开实验，选择 File → Export → Results（或导出为 .txt/.xls/.xlsx）后再上传。",
  sds: "这是 ABI SDS 系列仪器的实验工程文件（.sds），无法直接解析。请在 SDS 软件中 Export Results 为 .txt/.csv 后再上传。",
  pcrd: "这是 Bio-Rad CFX Maestro 的实验工程文件（.pcrd），无法直接解析。请在 CFX Maestro 中打开实验，使用 Export → Custom Export / Export All Data Sheets 导出为 .csv/.txt/.xlsx 后再上传。",
};

// --------------------------------------------------------------------------
// 基础小工具
// --------------------------------------------------------------------------
/** 列名/单元格归一化：去引号、压缩空白、转小写。 */
export function normCell(v: string | null | undefined): string {
  if (v === null || v === undefined) return "";
  let s = String(v).trim();
  // Python str.strip('"') / strip("'")：去掉两端所有引号字符
  s = s.replace(/^"+/, "").replace(/"+$/, "");
  s = s.replace(/^'+/, "").replace(/'+$/, "");
  s = s.trim();
  s = s.replace(/\s+/g, " ");
  return s.toLowerCase();
}

function hasCjk(s: string): boolean {
  for (const ch of s) {
    const code = ch.codePointAt(0) ?? 0;
    if (code >= 0x4e00 && code <= 0x9fff) return true;
  }
  return false;
}

/**
 * 文本解码：utf-8（含 BOM）→ gb18030（国产仪器 GBK）回退。
 * 对应 Python 版 _decode_text；浏览器 / Node（full-icu）均可用。
 */
export function decodeText(data: Uint8Array): string {
  try {
    return new TextDecoder("utf-8", { fatal: true }).decode(data);
  } catch {
    // fall through
  }
  try {
    return new TextDecoder("gb18030").decode(data);
  } catch {
    return new TextDecoder("utf-8").decode(data);
  }
}

// --------------------------------------------------------------------------
// 文本 → 二维网格
// --------------------------------------------------------------------------
/**
 * 单分隔符 CSV/TSV 解析，语义对齐 Python csv.reader
 * （引号仅在字段起始处生效，"" 为转义引号，支持引号内换行）。
 */
export function parseDelimited(text: string, sep: string): Grid {
  const rows: Grid = [];
  let row: string[] = [];
  let field = "";
  let inQuotes = false;
  let i = 0;
  const n = text.length;
  while (i < n) {
    const ch = text[i];
    if (inQuotes) {
      if (ch === '"') {
        if (i + 1 < n && text[i + 1] === '"') {
          field += '"';
          i += 2;
          continue;
        }
        inQuotes = false;
        i += 1;
        continue;
      }
      field += ch;
      i += 1;
      continue;
    }
    if (ch === '"' && field === "") {
      inQuotes = true;
      i += 1;
      continue;
    }
    if (ch === sep) {
      row.push(field);
      field = "";
      i += 1;
      continue;
    }
    if (ch === "\n") {
      row.push(field);
      rows.push(row);
      row = [];
      field = "";
      i += 1;
      continue;
    }
    field += ch;
    i += 1;
  }
  if (field !== "" || row.length > 0) {
    row.push(field);
    rows.push(row);
  }
  return rows;
}

/** 在候选分隔符中选出平均列数最多的一个（对齐 Python _choose_sep）。 */
export function chooseSep(lines: string[]): string {
  let bestSep = "\t";
  let bestCols = 0;
  const sampleLines = lines.slice(0, 50).filter((ln) => ln.trim() !== "").slice(0, 20);
  if (sampleLines.length === 0) return bestSep;
  for (const sep of ["\t", ",", ";"]) {
    const counts: number[] = [];
    for (const ln of sampleLines) {
      try {
        const parsed = parseDelimited(ln, sep);
        counts.push(parsed.length > 0 ? parsed[0].length : 1);
      } catch {
        counts.push(1);
      }
    }
    const avg = counts.reduce((a, b) => a + b, 0) / counts.length;
    if (avg > bestCols) {
      bestSep = sep;
      bestCols = avg;
    }
  }
  return bestSep;
}

/** 文本 → 网格：统一换行、嗅探分隔符、逐格 trim（对齐 Python _text_to_grid）。 */
export function textToGrid(text: string): Grid {
  const normalized = text.replace(/\r\n/g, "\n").replace(/\r/g, "\n");
  const lines = normalized.split("\n");
  const sep = chooseSep(lines);
  return parseDelimited(normalized, sep).map((row) => row.map((c) => c.trim()));
}

// --------------------------------------------------------------------------
// HTML 伪 xls → 网格（对齐 Python _html_to_grid；轻量实现，不依赖 DOM）
// --------------------------------------------------------------------------
function decodeHtmlEntities(s: string): string {
  return s
    .replace(/&#(\d+);/g, (_, d) => String.fromCodePoint(Number(d)))
    .replace(/&#x([0-9a-fA-F]+);/g, (_, h) => String.fromCodePoint(parseInt(h, 16)))
    .replace(/&nbsp;/g, " ")
    .replace(/&quot;/g, '"')
    .replace(/&#39;|&apos;/g, "'")
    .replace(/&lt;/g, "<")
    .replace(/&gt;/g, ">")
    .replace(/&amp;/g, "&");
}

/** 从 HTML 文本提取所有 <table> 为网格数组。 */
export function htmlTablesToGrids(text: string): Grid[] {
  const grids: Grid[] = [];
  const tableRe = /<table\b[\s\S]*?<\/table>/gi;
  let tm: RegExpExecArray | null;
  while ((tm = tableRe.exec(text)) !== null) {
    const tableHtml = tm[0];
    const grid: Grid = [];
    const trRe = /<tr\b[\s\S]*?<\/tr>/gi;
    let rm: RegExpExecArray | null;
    while ((rm = trRe.exec(tableHtml)) !== null) {
      const rowHtml = rm[0];
      const cells: string[] = [];
      const cellRe = /<t[dh]\b[^>]*>([\s\S]*?)<\/t[dh]>/gi;
      let cm: RegExpExecArray | null;
      while ((cm = cellRe.exec(rowHtml)) !== null) {
        const inner = cm[1].replace(/<[^>]+>/g, "");
        cells.push(decodeHtmlEntities(inner).trim());
      }
      if (cells.length > 0) grid.push(cells);
    }
    if (grid.length > 0) grids.push(grid);
  }
  return grids;
}

/** 判断文本是否是 HTML 伪 xls（对齐 Python _sniff_file_type 的 html 分支）。 */
export function looksLikeHtml(text: string): boolean {
  const head = text.slice(0, 1024).replace(/^\s+/, "").toLowerCase();
  return (
    head.startsWith("<") &&
    (head.includes("<table") ||
      head.includes("<html") ||
      head.includes("<!doctype html") ||
      head.includes("<meta"))
  );
}

/** HTML 伪 xls：取含结果表头得分最高的表（对齐 Python _html_to_grid）。 */
export function htmlToGrid(name: string, text: string): Grid {
  const tables = htmlTablesToGrids(text);
  let best: { score: number; grid: Grid } | null = null;
  for (const grid of tables) {
    for (let i = 0; i < Math.min(40, grid.length); i++) {
      const { score } = scoreRow(grid[i]);
      if (best === null || score > best.score) {
        best = { score, grid };
      }
    }
  }
  if (best === null || best.score < 2) {
    throw new InstrumentParseError(
      `${name} 是 HTML 格式的结果文件，但未找到 qPCR 结果表。` +
        "请改用仪器软件导出为 .txt/.csv/.xlsx，或使用手动解析。"
    );
  }
  return best.grid;
}

// --------------------------------------------------------------------------
// 表头行扫描与仪器识别
// --------------------------------------------------------------------------
/** 给一行打分：归一化后精确命中同义词表的字段数。 */
export function scoreRow(cells: string[]): { score: number; fields: Set<Field> } {
  const normCells = cells.map((c) => normCell(c));
  const fields = new Set<Field>();
  for (const field of Object.keys(SYNONYMS) as Field[]) {
    const syns = SYNONYMS[field];
    for (const c of normCells) {
      if (c !== "" && syns.includes(c)) {
        fields.add(field);
        break;
      }
    }
  }
  return { score: fields.size, fields };
}

/** 扫描前 80 行，找同义词命中最多、且含 ct + (well 或 sample) 的表头行。 */
export function detectHeaderRow(grid: Grid): { headerIdx: number; fields: Set<Field> } {
  let best: { score: number; idx: number; fields: Set<Field> } | null = null;
  for (let i = 0; i < Math.min(80, grid.length); i++) {
    const { score, fields } = scoreRow(grid[i]);
    if (!fields.has("ct")) continue;
    if (!fields.has("well") && !fields.has("sample")) continue;
    if (best === null || score > best.score) {
      best = { score, idx: i, fields };
    }
  }
  if (best === null) {
    throw new InstrumentParseError(
      "未能定位结果表头行（应同时包含 Ct/Cq/Cp 与 孔位或样本 列）。" +
        "请确认文件是仪器导出的结果表；若表头格式特殊，请使用下方手动解析。"
    );
  }
  return { headerIdx: best.idx, fields: best.fields };
}

/** 仪器规则链分类（顺序与 Python 版一致）。 */
export function classifyInstrument(headerCells: string[]): string {
  const s = new Set(headerCells.map((c) => normCell(c)));
  s.delete("");
  if (s.has("cp") && (s.has("include") || s.has("color") || s.has("status"))) {
    return "Roche LightCycler 480";
  }
  if (s.has("content") || s.has("fluor") || s.has("fluorophore")) {
    return "Bio-Rad CFX";
  }
  if (s.has("task") || s.has("well position") || s.has("target name")) {
    return "ABI/Thermo QuantStudio";
  }
  if (s.has("gene")) {
    return "Roche LightCycler 96";
  }
  if ((s.has("no.") || s.has("no")) && s.has("ct")) {
    return "Qiagen Rotor-Gene Q";
  }
  if (s.has("dye")) {
    return "Agilent AriaMx";
  }
  if (s.has("well") && s.has("target") && s.has("ct")) {
    return "ABI/Thermo QuantStudio (v2.x)";
  }
  for (const c of s) {
    if (hasCjk(c)) return "国产仪器/通用格式（中文表头）";
  }
  return "通用格式（按列名同义词识别）";
}

/** 列名 → well/sample/target/ct 映射。ct 优先精确同义词，兜底均值列。 */
export function mapColumns(columns: string[]): Partial<Record<Field, string>> {
  const normMap = columns.map((c) => ({ col: c, norm: normCell(c) }));
  const mapping: Partial<Record<Field, string>> = {};
  const used = new Set<string>();
  for (const field of ["ct", "target", "sample", "well"] as Field[]) {
    for (const syn of SYNONYMS[field]) {
      const hit = normMap.find((e) => e.norm === syn && !used.has(e.col));
      if (hit !== undefined) {
        mapping[field] = hit.col;
        used.add(hit.col);
        break;
      }
    }
  }
  if (mapping.ct === undefined) {
    // 兜底：只有 Cq Mean 等列
    for (const syn of CT_FALLBACK) {
      const hit = normMap.find((e) => e.norm === syn && !used.has(e.col));
      if (hit !== undefined) {
        mapping.ct = hit.col;
        used.add(hit.col);
        break;
      }
    }
  }
  return mapping;
}

/** 在 raw 列中按归一化名称找列（对齐 Python _get_col，后者覆盖前者）。 */
function getCol(raw: RawTable, ...names: string[]): string | null {
  const norm = new Map<string, string>();
  for (const c of raw.columns) norm.set(normCell(c), c);
  for (const n of names) {
    const hit = norm.get(normCell(n));
    if (hit !== undefined) return hit;
  }
  return null;
}

/** 统一 Ct 清洗：未检出 → null；欧区小数逗号 → 小数点。 */
export function cleanCt(v: string): number | null {
  const s = normCell(v);
  if (UNDETERMINED.has(s)) return null;
  let raw = String(v).trim();
  if (/^-?\d{1,3},\d+$/.test(raw)) {
    // 23,5729 → 23.5729
    raw = raw.replace(",", ".");
  }
  const num = Number(raw);
  return Number.isNaN(num) ? null : num;
}

/** 表头去重：空名补 col_j，重名加 _N 后缀（对齐 Python _dedup_header）。 */
export function dedupHeader(header: string[]): string[] {
  const out: string[] = [];
  const seen = new Map<string, number>();
  for (let j = 0; j < header.length; j++) {
    let h = String(header[j]).trim();
    if (h === "") h = `col_${j}`;
    if (seen.has(h)) {
      const next = (seen.get(h) ?? 0) + 1;
      seen.set(h, next);
      h = `${h}_${next}`;
    } else {
      seen.set(h, 0);
    }
    out.push(h);
  }
  return out;
}

// --------------------------------------------------------------------------
// 网格 → 标准化长表
// --------------------------------------------------------------------------
/**
 * 网格 + 表头行号 → 原始表：
 * 跳过全空行、截断/补齐到表头宽度、去掉自动命名且全空的列
 * （对齐 Python _grid_to_raw_df）。
 */
export function gridToRawTable(grid: Grid, headerIdx: number): RawTable {
  const header = dedupHeader(grid[headerIdx]);
  let rows: string[][] = [];
  for (const row of grid.slice(headerIdx + 1)) {
    if (!row.some((c) => String(c).trim() !== "")) continue;
    const padded = row.slice(0, header.length);
    while (padded.length < header.length) padded.push("");
    rows.push(padded);
  }
  // 去掉全空列（Bio-Rad 首列前多余分隔符产生的空列）
  const keepIdx: number[] = [];
  for (let j = 0; j < header.length; j++) {
    const allEmpty = rows.every((r) => String(r[j]).trim() === "");
    if (!(allEmpty && header[j].startsWith("col_"))) keepIdx.push(j);
  }
  return {
    columns: keepIdx.map((j) => header[j]),
    rows: rows.map((r) => keepIdx.map((j) => r[j])),
  };
}

/** LC480 元信息行里的 Filter Combination 通道标识，如 465-510。 */
export function lc480Channel(grid: Grid, headerIdx: number): string | null {
  for (const row of grid.slice(0, headerIdx)) {
    const line = row.map((c) => String(c)).join(" ");
    const m = /filter combination[:\s]*([0-9]{3}\s*[-–]\s*[0-9]{3})/i.exec(line);
    if (m) return m[1].replace(/\s+/g, "");
  }
  return null;
}

export interface StandardizeResult {
  rows: LongRow[];
  needsTargetName: boolean;
}

/** 原始表 → 标准化长表（对齐 Python _standardize，含仪器特定过滤）。 */
export function standardize(
  raw: RawTable,
  instrument: string,
  sourceFile: string,
  defaultTarget: string | null,
  manualMap?: Partial<Record<Field, string>>
): StandardizeResult {
  const mapping = manualMap ?? mapColumns(raw.columns);

  if (mapping.ct === undefined || !raw.columns.includes(mapping.ct)) {
    throw new InstrumentParseError(
      "未找到 Ct/Cq/Cp 列。请确认导出的结果表包含 Ct 值，或使用手动解析指定列。"
    );
  }
  if (mapping.sample === undefined && mapping.well === undefined) {
    throw new InstrumentParseError(
      "未找到 样本/孔位 列。请确认导出的结果表包含样本名或孔位，或使用手动解析。"
    );
  }

  const colIdx = (name: string): number => raw.columns.indexOf(name);
  const cellAt = (row: string[], name: string | undefined): string =>
    name === undefined ? "" : row[colIdx(name)] ?? "";

  const needsTarget = mapping.target === undefined;
  const fallbackTarget = defaultTarget || `待命名靶标(${sourceFile})`;

  // 仪器特定过滤所需的列（按原始列名解析一次）
  const taskCol = getCol(raw, "Task");
  const contentCol = getCol(raw, "Content");
  const omitCol = getCol(raw, "Omit");
  const typeCol = getCol(raw, "Type");

  const out: LongRow[] = [];
  for (const row of raw.rows) {
    // ---- 仪器特定过滤（顺序与 Python 版一致）----
    if (taskCol !== null && normCell(cellAt(row, taskCol)) === "ntc") continue;
    if (
      contentCol !== null &&
      ["ntc", "no template control"].includes(normCell(cellAt(row, contentCol)))
    ) {
      continue;
    }
    if (
      omitCol !== null &&
      ["true", "1", "yes"].includes(normCell(cellAt(row, omitCol)))
    ) {
      continue;
    }
    if (
      typeCol !== null &&
      instrument === "Qiagen Rotor-Gene Q" &&
      ["ntc", "no template control"].includes(normCell(cellAt(row, typeCol)))
    ) {
      continue;
    }

    const well = mapping.well !== undefined ? cellAt(row, mapping.well).trim() : "";
    const sample =
      mapping.sample !== undefined ? cellAt(row, mapping.sample).trim() : well;
    const target =
      mapping.target !== undefined ? cellAt(row, mapping.target).trim() : fallbackTarget;
    const ct = cleanCt(cellAt(row, mapping.ct));

    // 空样本行（未加样的空孔）
    if (["", "none", "nan"].includes(normCell(sample))) continue;

    out.push({ well, sample, target, ct, source_file: sourceFile });
  }
  return { rows: out, needsTargetName: needsTarget };
}

/** 自动模式：网格 → (长表, 报告)（对齐 Python _parse_grid）。 */
export function parseGrid(grid: Grid, sourceFile: string): ParseResult {
  const { headerIdx } = detectHeaderRow(grid);
  const instrument = classifyInstrument(grid[headerIdx]);
  const raw = gridToRawTable(grid, headerIdx);
  if (raw.rows.length === 0) {
    throw new InstrumentParseError("表头之后没有数据行，请确认导出了结果数据。");
  }

  // 单通道文件的默认 target 命名
  let defaultTarget: string | null = null;
  if (instrument === "Roche LightCycler 480") {
    const ch = lc480Channel(grid, headerIdx);
    defaultTarget = ch ? `LC480通道${ch}` : null;
  } else if (instrument === "Qiagen Rotor-Gene Q") {
    defaultTarget = sourceFile.includes(".")
      ? sourceFile.slice(0, sourceFile.lastIndexOf("."))
      : sourceFile;
  }

  const { rows, needsTargetName } = standardize(raw, instrument, sourceFile, defaultTarget);
  const targets = [...new Set(rows.map((r) => r.target))].sort();
  const report: FileReport = {
    file: sourceFile,
    instrument,
    status: "ok",
    rows: rows.length,
    targets,
    needs_target_name: needsTargetName,
    default_target: defaultTarget || "",
    message: `识别为 ${instrument}，解析 ${rows.length} 行。`,
  };
  return { rows, report };
}

// --------------------------------------------------------------------------
// 文本文件入口（txt / csv / HTML 伪 xls）
// --------------------------------------------------------------------------
/**
 * 解析单个文本类文件（txt/csv，以及 HTML 伪 xls）。
 * 调用方负责把字节转码为字符串（可用导出的 decodeText 做 utf-8→gb18030 嗅探）；
 * 本函数内做：二进制工程文件扩展名拦截、HTML 伪 xls 表格提取、
 * 分隔符嗅探、表头定位与仪器识别。
 */
export function parseTextFile(filename: string, text: string): ParseResult {
  const lower = filename.toLowerCase();
  const ext = lower.includes(".") ? lower.slice(lower.lastIndexOf(".") + 1) : "";
  if (ext in BINARY_HINTS) {
    throw new InstrumentParseError(BINARY_HINTS[ext]);
  }
  const grid = looksLikeHtml(text) ? htmlToGrid(filename, text) : textToGrid(text);
  return parseGrid(grid, filename);
}

/**
 * 多文件自动解析入口（对齐 Python parse_instrument_files）。
 * 解析失败的文件不会中断其他文件，错误写进 reports。
 */
export function parseTextFiles(files: Array<{ name: string; text: string }>): MultiParseResult {
  const rows: LongRow[] = [];
  const reports: FileReport[] = [];
  for (const { name, text } of files) {
    try {
      const { rows: fileRows, report } = parseTextFile(name, text);
      rows.push(...fileRows);
      reports.push(report);
    } catch (e) {
      if (e instanceof InstrumentParseError) {
        reports.push({
          file: name,
          instrument: "未知",
          status: "error",
          rows: 0,
          targets: [],
          needs_target_name: false,
          default_target: "",
          message: e.message,
        });
      } else {
        reports.push({
          file: name,
          instrument: "未知",
          status: "error",
          rows: 0,
          targets: [],
          needs_target_name: false,
          default_target: "",
          message: `解析时出现未预期错误：${String(e)}。可尝试手动解析。`,
        });
      }
    }
  }
  return { rows, reports };
}

// --------------------------------------------------------------------------
// 手动模式
// --------------------------------------------------------------------------
/**
 * 手动模式：用户指定表头行（0 基）与各列名（对齐 Python parse_manual，
 * 但直接接收网格——文本/Excel 的网格化在调用方完成）。
 */
export function parseManualGrid(
  grid: Grid,
  sourceFile: string,
  headerRow: number,
  cols: ManualColumns,
  defaultTarget?: string
): ParseResult {
  if (!(0 <= headerRow && headerRow < grid.length)) {
    throw new InstrumentParseError(
      `表头行号 ${headerRow + 1} 超出范围（共 ${grid.length} 行）。`
    );
  }
  const raw = gridToRawTable(grid, headerRow);
  const manualMap: Partial<Record<Field, string>> = {};
  const candidates: Array<[Field, string | undefined]> = [
    ["well", cols.well],
    ["sample", cols.sample],
    ["target", cols.target],
    ["ct", cols.ct],
  ];
  for (const [field, col] of candidates) {
    if (col !== undefined && col !== "" && raw.columns.includes(col)) {
      manualMap[field] = col;
    }
  }
  const instrument = classifyInstrument(grid[headerRow]);
  const { rows, needsTargetName } = standardize(
    raw,
    instrument,
    sourceFile,
    defaultTarget ?? null,
    manualMap
  );
  const targets = [...new Set(rows.map((r) => r.target))].sort();
  const report: FileReport = {
    file: sourceFile,
    instrument: `${instrument}（手动指定）`,
    status: "ok",
    rows: rows.length,
    targets,
    needs_target_name: needsTargetName,
    default_target: defaultTarget || "",
    message: `手动解析成功，共 ${rows.length} 行。`,
  };
  return { rows, report };
}

// --------------------------------------------------------------------------
// 下游辅助：分组与宽表
// --------------------------------------------------------------------------
/**
 * 样本 → 分组 默认映射。
 * rule:
 *   identity        每个样本独立成组
 *   strip_suffix    前缀分组：去掉尾部重复编号（如 _1、-2、空格3、.1）
 *   prefix_sep      前缀分组：取第一个 _ - 空格 之前的部分
 */
export function autoGroupSamples(
  samples: string[],
  rule: "identity" | "strip_suffix" | "prefix_sep" = "identity"
): Record<string, string> {
  const mapping: Record<string, string> = {};
  for (const raw of samples) {
    const s = String(raw);
    let g: string;
    if (rule === "strip_suffix") {
      g = s.replace(/[\s_\-\.]*\d+$/, "").trim() || s;
    } else if (rule === "prefix_sep") {
      g = s.split(/[_\-\s]/, 2)[0].trim() || s;
    } else {
      g = s;
    }
    mapping[s] = g;
  }
  return mapping;
}

/** 出现在所有样本中的 target（候选管家基因）置顶，返回排序后的 target 列表。 */
export function targetsInAllSamples(rows: LongRow[]): string[] {
  if (rows.length === 0) return [];
  const samples = new Set(rows.map((r) => String(r.sample)));
  const allTargets = [...new Set(rows.map((r) => String(r.target)))];
  const coverage = (t: string): number =>
    new Set(rows.filter((r) => String(r.target) === t).map((r) => String(r.sample))).size;
  return allTargets.sort((a, b) => {
    const ca = coverage(a);
    const cb = coverage(b);
    const aAll = ca === samples.size ? 0 : 1;
    const bAll = cb === samples.size ? 0 : 1;
    if (aAll !== bAll) return aAll - bAll;
    if (ca !== cb) return cb - ca;
    return a < b ? -1 : a > b ? 1 : 0;
  });
}

export interface WideTable {
  /** 列名：sample + 各 target（target 列按名称排序，对齐 pandas unstack）。 */
  columns: string[];
  /** 每行一个样本；ct 为同一样本同一靶标多孔均值，全未检出为 null。 */
  rows: Array<Record<string, string | number | null>>;
}

/** 长表 → 宽表：每行一个样本，每个 target 一列（同样本同靶标多孔取均值）。 */
export function pivotToWide(rows: LongRow[]): WideTable {
  const sampleSet = new Set<string>();
  const targetSet = new Set<string>();
  const sums = new Map<string, { sum: number; n: number }>();
  for (const r of rows) {
    const s = String(r.sample);
    const t = String(r.target);
    sampleSet.add(s);
    targetSet.add(t);
    if (r.ct !== null) {
      const key = `${s}\u0000${t}`;
      const acc = sums.get(key) ?? { sum: 0, n: 0 };
      acc.sum += r.ct;
      acc.n += 1;
      sums.set(key, acc);
    }
  }
  // pandas groupby/unstack 默认按 key 排序
  const samples = [...sampleSet].sort();
  const targets = [...targetSet].sort();
  const wideRows = samples.map((s) => {
    const rec: Record<string, string | number | null> = { sample: s };
    for (const t of targets) {
      const acc = sums.get(`${s}\u0000${t}`);
      rec[t] = acc !== undefined && acc.n > 0 ? acc.sum / acc.n : null;
    }
    return rec;
  });
  return { columns: ["sample", ...targets], rows: wideRows };
}
