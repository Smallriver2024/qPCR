/**
 * templates.ts — 新手友好的下载模板
 *
 * - downloadOrganizedTemplate：「上传整理表」用的 xlsx 模板
 *   （第一个 sheet 为数据模板，含示例数据；第二个 sheet 为填写说明。
 *    注意 readTableFile 只读第一个工作表，数据模板必须放最前。）
 * - downloadInstrumentExample：「仪器原始文件」用的示例导出文件
 *   （QuantStudio v1.x 风格 txt / Bio-Rad CFX 风格 csv，列结构对齐
 *    parsers.ts 的 classifyInstrument + SYNONYMS，可直接上传体验完整流程。）
 *
 * 示例数据与 Home.tsx 的 PASTE_PLACEHOLDER 同源：NC 3 重复 + TreatmentA 3 重复，
 * 管家基因 GAPDH，目标基因 MYC，计算后 TreatmentA 相对表达约 2.46。
 */
import * as XLSX from "xlsx";
import { downloadBlob } from "./exporters.ts";

// --------------------------------------------------------------------------
// 共享示例数据（sample → group 映射用 _1 后缀规则即可自动分组）
// --------------------------------------------------------------------------
const SAMPLES: Array<{ sample: string; group: string; gapdh: number; myc: number }> = [
  { sample: "NC_1", group: "NC", gapdh: 12.1, myc: 23.2 },
  { sample: "NC_2", group: "NC", gapdh: 12.3, myc: 22.8 },
  { sample: "NC_3", group: "NC", gapdh: 12.2, myc: 23.0 },
  { sample: "TreatmentA_1", group: "TreatmentA", gapdh: 12.4, myc: 21.9 },
  { sample: "TreatmentA_2", group: "TreatmentA", gapdh: 12.5, myc: 21.6 },
  { sample: "TreatmentA_3", group: "TreatmentA", gapdh: 12.3, myc: 21.7 },
];

// --------------------------------------------------------------------------
// 整理表模板（xlsx：数据模板 + 填写说明）
// --------------------------------------------------------------------------
export function downloadOrganizedTemplate(): void {
  const wb = XLSX.utils.book_new();

  // Sheet 1：数据模板（readTableFile 只读第一个工作表）
  const dataAoa: Array<Array<string | number>> = [
    ["group", "sample", "GAPDH", "MYC"],
    ...SAMPLES.map((s) => [s.group, s.sample, s.gapdh, s.myc] as Array<string | number>),
  ];
  const wsData = XLSX.utils.aoa_to_sheet(dataAoa);
  wsData["!cols"] = [{ wch: 14 }, { wch: 14 }, { wch: 10 }, { wch: 10 }];
  XLSX.utils.book_append_sheet(wb, wsData, "数据模板");

  // Sheet 2：填写说明
  const guideAoa = [
    ["qPCR ΔΔCt 计算器 —— 整理表填写说明"],
    [""],
    ["1. 第一行必须是列名，请勿删除或改名前两列（group、sample）。"],
    ["2. group 列：填分组名（如 NC、TreatmentA），同一组的重复样本填相同的组名。"],
    ["3. sample 列：填样本名，每个样本一行；该列可留空或整列删除。"],
    ["4. 第三列起为各基因的 Ct 值，列名即基因名（如 GAPDH、MYC），可任意增删基因列。"],
    ["5. 至少需要 2 个基因列：一个目标基因 + 一个管家基因（如 GAPDH、ACTB、18S）。"],
    ["6. Ct 未检出的孔请留空，不要填 0 或文字。"],
    ["7. 模板中的 6 行是示例数据（NC / TreatmentA 各 3 个重复），可直接上传试算，之后替换为自己的数据。"],
    ["8. 填好后保存本文件，到网站「上传整理表」中上传即可；也支持另存为 .csv 上传。"],
  ];
  const wsGuide = XLSX.utils.aoa_to_sheet(guideAoa);
  wsGuide["!cols"] = [{ wch: 88 }];
  XLSX.utils.book_append_sheet(wb, wsGuide, "填写说明");

  const out = XLSX.write(wb, { type: "array", bookType: "xlsx" }) as ArrayBuffer;
  downloadBlob(
    "qPCR整理表模板.xlsx",
    new Blob([out], {
      type: "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
    }),
  );
}

// --------------------------------------------------------------------------
// 仪器示例文件
// --------------------------------------------------------------------------
/** ABI/Thermo QuantStudio v1.x 风格结果导出（制表符分隔 txt，含元信息头）。 */
function quantstudioExample(): string {
  const header = [
    "Well", "Well Position", "Sample Name", "Target Name",
    "Task", "Reporter", "Quencher", "CT", "Ct Mean",
  ];
  const lines: string[] = [
    "* Block Type = 96-Well Block (0.2 mL)",
    "* Chemistry = TAQMAN",
    "* Experiment File Name = demo_experiment.eds",
    "* Experiment Run End Time = 2026-07-01 10:30:00",
    "* Instrument Type = QuantStudio(TM) 5",
    "",
    "[Results]",
    header.join("\t"),
  ];
  let well = 0;
  const pos = (i: number): string =>
    `${String.fromCharCode(65 + Math.floor(i / 12))}${(i % 12) + 1}`;
  for (const s of SAMPLES) {
    for (const [target, ct] of [["GAPDH", s.gapdh], ["MYC", s.myc]] as Array<[string, number]>) {
      well += 1;
      lines.push(
        [String(well), pos(well - 1), s.sample, target, "UNKNOWN", "FAM", "NFQ", String(ct), String(ct)].join("\t"),
      );
    }
  }
  return lines.join("\r\n");
}

/** Bio-Rad CFX Maestro 风格结果导出（逗号分隔 csv）。 */
function bioradExample(): string {
  const header = ["Well", "Fluor", "Content", "Sample", "Target", "Cq", "Starting Quantity (SQ)"];
  const lines = [header.join(",")];
  let well = 0;
  const pos = (i: number): string =>
    `${String.fromCharCode(65 + Math.floor(i / 12))}${String((i % 12) + 1).padStart(2, "0")}`;
  for (const s of SAMPLES) {
    for (const [target, ct] of [["GAPDH", s.gapdh], ["MYC", s.myc]] as Array<[string, number]>) {
      well += 1;
      lines.push(
        [pos(well - 1), "SYBR", `Unkn-${String(well).padStart(2, "0")}`, s.sample, target, String(ct), ""].join(","),
      );
    }
  }
  return lines.join("\r\n");
}

export type InstrumentExampleKind = "quantstudio" | "biorad";

/** 下载仪器示例文件（可直接上传到「仪器原始文件」体验完整流程）。 */
export function downloadInstrumentExample(kind: InstrumentExampleKind): void {
  if (kind === "quantstudio") {
    downloadBlob(
      "示例_QuantStudio结果导出.txt",
      new Blob([quantstudioExample()], { type: "text/plain;charset=utf-8" }),
    );
  } else {
    downloadBlob(
      "示例_BioRad_CFX结果导出.csv",
      new Blob(["\uFEFF" + bioradExample()], { type: "text/csv;charset=utf-8" }),
    );
  }
}
