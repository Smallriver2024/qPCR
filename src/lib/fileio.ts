/**
 * fileio.ts — 浏览器端文件读取：字节 → 二维网格 / 整理表
 *
 * 设计简报 §11：
 *   Excel 读取用 SheetJS 读 ArrayBuffer → 多 sheet 用 parsers 的 scoreRow
 *   打分选表 → 二维网格喂 parseGrid；老 .xls 同样走 SheetJS；
 *   GBK 文本用 decodeText。
 */
import * as XLSX from "xlsx";
import {
  decodeText,
  dedupHeader,
  htmlToGrid,
  InstrumentParseError,
  looksLikeHtml,
  parseGrid,
  parseTextFile,
  scoreRow,
  textToGrid,
  BINARY_HINTS,
  type Grid,
  type ParseResult,
} from "./parsers.ts";
import type { PreparedTable } from "./analysis.ts";

const EXCEL_EXTS = new Set(["xls", "xlsx", "xlsm"]);

export function extOf(name: string): string {
  const lower = name.toLowerCase();
  return lower.includes(".") ? lower.slice(lower.lastIndexOf(".") + 1) : "";
}

/** SheetJS worksheet → 字符串网格。 */
function sheetToGrid(ws: XLSX.WorkSheet): Grid {
  const aoa = XLSX.utils.sheet_to_json<unknown[]>(ws, {
    header: 1,
    raw: false,
    defval: "",
  });
  return aoa.map((row) => (Array.isArray(row) ? row.map((c) => String(c ?? "")) : []));
}

/**
 * Excel 字节 → 网格：对每个 sheet 的前 40 行用 scoreRow 打分，
 * 选得分最高的工作表（对齐旧版 parsers.load_grid 的多 sheet 逻辑）。
 */
export function excelToGrid(name: string, buf: ArrayBuffer): Grid {
  let wb: XLSX.WorkBook;
  try {
    wb = XLSX.read(buf, { type: "array" });
  } catch {
    throw new InstrumentParseError(
      `${name} 无法作为 Excel 读取。若这是早期仪器软件导出的伪 xls，请改用 txt/csv 导出，或使用手动解析。`,
    );
  }
  let best: { score: number; grid: Grid } | null = null;
  for (const sheetName of wb.SheetNames) {
    const grid = sheetToGrid(wb.Sheets[sheetName]);
    for (let i = 0; i < Math.min(40, grid.length); i++) {
      const { score } = scoreRow(grid[i]);
      if (best === null || score > best.score) best = { score, grid };
    }
  }
  if (best === null || best.score < 2) {
    throw new InstrumentParseError(
      `${name} 的各工作表中均未找到 qPCR 结果表头（应包含 Ct/Cq/Cp 与 孔位/样本 列）。可尝试手动解析。`,
    );
  }
  return best.grid;
}

/** 任意仪器文件字节 → 网格（手动解析用；二进制工程文件给出中文导出指引）。 */
export function fileToGrid(name: string, buf: ArrayBuffer): Grid {
  const ext = extOf(name);
  if (ext in BINARY_HINTS) throw new InstrumentParseError(BINARY_HINTS[ext]);
  if (EXCEL_EXTS.has(ext)) return excelToGrid(name, buf);
  const text = decodeText(new Uint8Array(buf));
  return looksLikeHtml(text) ? htmlToGrid(name, text) : textToGrid(text);
}

/** 自动模式：解析单个仪器原始文件（txt/csv/tsv/xls/xlsx/eds/sds/pcrd）。 */
export function parseInstrumentFile(name: string, buf: ArrayBuffer): ParseResult {
  const ext = extOf(name);
  if (EXCEL_EXTS.has(ext)) return parseGrid(excelToGrid(name, buf), name);
  const text = ext in BINARY_HINTS ? "" : decodeText(new Uint8Array(buf));
  return parseTextFile(name, text); // 二进制扩展名在 parseTextFile 内拦截并给出指引
}

// --------------------------------------------------------------------------
// 整理好的表格（复制粘贴 / 上传整理表）→ PreparedTable
// --------------------------------------------------------------------------
/** 网格 → 整理表：首行为列名（去重），跳过全空行。 */
function gridToPrepared(grid: Grid): PreparedTable {
  const nonempty = grid.filter((row) => row.some((c) => String(c).trim() !== ""));
  if (nonempty.length < 2) {
    throw new InstrumentParseError("未读到有效数据：至少需要 1 行列名 + 1 行数据。");
  }
  const width = Math.max(...nonempty.map((r) => r.length));
  const header = dedupHeader(
    Array.from({ length: width }, (_, j) => String(nonempty[0][j] ?? "")),
  );
  const rows = nonempty.slice(1).map((row) => {
    const rec: Record<string, string | number | null> = {};
    for (let j = 0; j < header.length; j++) rec[header[j]] = String(row[j] ?? "").trim();
    return rec;
  });
  return { columns: header, rows };
}

/** 粘贴文本 → 整理表（分隔符嗅探复用 parsers.textToGrid）。 */
export function readPastedTable(text: string): PreparedTable {
  if (!text.trim()) {
    throw new InstrumentParseError("请先粘贴数据。至少包含：group 列 + 2 列 Ct 数值（目标与管家）。");
  }
  return gridToPrepared(textToGrid(text));
}

/** 上传整理表（csv/tsv/txt/xls/xlsx）→ 整理表：取第一个工作表，首行为列名。 */
export function readTableFile(name: string, buf: ArrayBuffer): PreparedTable {
  const ext = extOf(name);
  if (EXCEL_EXTS.has(ext)) {
    let wb: XLSX.WorkBook;
    try {
      wb = XLSX.read(buf, { type: "array" });
    } catch {
      throw new InstrumentParseError(`${name} 无法作为 Excel 读取，请改存为 .xlsx 或 .csv 后重试。`);
    }
    const ws = wb.Sheets[wb.SheetNames[0]];
    if (!ws) throw new InstrumentParseError(`${name} 中没有工作表。`);
    return gridToPrepared(sheetToGrid(ws));
  }
  return gridToPrepared(textToGrid(decodeText(new Uint8Array(buf))));
}
