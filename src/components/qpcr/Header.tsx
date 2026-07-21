/**
 * Header.tsx — 固定顶栏（design/kimi-ui-brief.md §4）+ 使用指南 Dialog
 */
import { useCallback, useEffect, useState, type ReactNode } from "react";
import { KimiButton } from "./ui.tsx";

function GuideDialog({ open, onClose }: { open: boolean; onClose: () => void }) {
  useEffect(() => {
    if (!open) return;
    const onKey = (e: KeyboardEvent) => {
      if (e.key === "Escape") onClose();
    };
    window.addEventListener("keydown", onKey);
    return () => window.removeEventListener("keydown", onKey);
  }, [open, onClose]);

  if (!open) return null;
  return (
    <div
      className="fixed inset-0 z-50 flex items-center justify-center bg-black/30 p-4"
      onClick={onClose}
      role="presentation"
    >
      <div
        role="dialog"
        aria-modal="true"
        aria-label="使用指南"
        className="max-h-[80vh] w-full max-w-[560px] overflow-auto rounded-lg border border-sep-s1 bg-white p-5"
        onClick={(e) => e.stopPropagation()}
      >
        <div className="mb-4 flex items-center justify-between">
          <h2 className="text-base font-medium leading-6 text-label-primary">使用指南</h2>
          <KimiButton variant="outline" onClick={onClose} aria-label="关闭">
            关闭
          </KimiButton>
        </div>
        <GuideContent />
      </div>
    </div>
  );
}

function GuideSection({ title, children }: { title: string; children: ReactNode }) {
  return (
    <div className="mb-4">
      <h3 className="mb-1.5 text-sm font-medium leading-5 text-label-primary">{title}</h3>
      <div className="space-y-1 text-sm leading-5 text-label-secondary">{children}</div>
    </div>
  );
}

function GuideContent() {
  return (
    <div>
      <GuideSection title="1. 准备数据（三选一）">
        <p>
          <strong>复制粘贴</strong>：粘贴整理好的表格，首行为列名，需包含 group 列 +
          至少两列 Ct 数值（目标基因与管家基因），支持制表符/逗号分隔。
        </p>
        <p>
          <strong>上传整理表</strong>：上传同样结构的 CSV / Excel 文件。
        </p>
        <p>
          <strong>仪器原始文件</strong>：直接上传仪器导出的原始结果文件，支持
          ABI/Thermo QuantStudio、Bio-Rad CFX、Roche LightCycler 480/96、Agilent
          AriaMx、Qiagen Rotor-Gene Q、Eppendorf 及国产仪器（天隆/博日/宏石）。
          分通道导出的仪器（如 LC480、Rotor-Gene）请一次上传全部通道文件。
        </p>
        <p>
          .eds / .sds / .pcrd 为二进制工程文件，请先在仪器软件中 Export 为
          txt / csv / xlsx 后再上传。
        </p>
      </GuideSection>
      <GuideSection title="2. 分析设置">
        <p>选择目标基因列、管家基因列，并填写对照组名称（与 group 列中的取值一致）。</p>
        <p>可选择 Welch t 检验或 Mann–Whitney U 检验，以及多重比较校正方式。</p>
      </GuideSection>
      <GuideSection title="3. 计算方法（ΔΔCt）">
        <p>ΔCt = 目标基因 Ct − 管家基因 Ct</p>
        <p>ΔΔCt = ΔCt − 对照组 ΔCt 均值</p>
        <p>相对表达量 = 2^(−ΔΔCt)</p>
      </GuideSection>
      <GuideSection title="4. 结果与导出">
        <p>
          计算后可查看逐样本结果、分组均值±SD、带显著性标注的柱状图与两两比较表，
          并导出 Excel（4 个 sheet）、两两比较 CSV 与柱状图 PNG。
        </p>
        <p>全部计算均在浏览器本地完成，数据不会上传到任何服务器。</p>
      </GuideSection>
    </div>
  );
}

export function Header() {
  const [guideOpen, setGuideOpen] = useState(false);
  const closeGuide = useCallback(() => setGuideOpen(false), []);
  return (
    <>
      <header className="fixed inset-x-0 top-0 z-40 h-16 border-b border-sep-s1 bg-white/70 backdrop-blur-[6px]">
        <div className="mx-auto flex h-full max-w-[880px] items-center justify-between px-4">
          <div className="flex items-center gap-3">
            <div className="flex h-8 w-8 items-center justify-center rounded-full bg-label-primary">
              <span className="text-sm font-semibold text-white">qP</span>
            </div>
            <div>
              <div className="text-base font-medium leading-6 text-label-primary">
                qPCR ΔΔCt 计算器
              </div>
              <div className="text-xs leading-[18px] text-label-tertiary">相对表达量分析</div>
            </div>
          </div>
          <KimiButton variant="outline" onClick={() => setGuideOpen(true)}>
            使用指南
          </KimiButton>
        </div>
      </header>
      <GuideDialog open={guideOpen} onClose={closeGuide} />
    </>
  );
}
