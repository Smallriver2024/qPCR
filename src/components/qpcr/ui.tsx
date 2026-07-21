/**
 * ui.tsx — Kimi 风格基础组件（design/kimi-ui-brief.md §3/§5/§6/§7）
 */
import { forwardRef, type ButtonHTMLAttributes, type ReactNode, type SelectHTMLAttributes } from "react";

// --------------------------------------------------------------------------
// 按钮（§5）：primary 黑底白字 / secondary f1 底 / outline 透明+s1 边
// --------------------------------------------------------------------------
type ButtonVariant = "primary" | "secondary" | "outline" | "dangerText";
type ButtonSize = "default" | "hero";

interface KimiButtonProps extends ButtonHTMLAttributes<HTMLButtonElement> {
  variant?: ButtonVariant;
  size?: ButtonSize;
  block?: boolean;
}

const variantCls: Record<ButtonVariant, string> = {
  primary:
    "bg-label-primary text-white hover:opacity-85 disabled:text-label-quaternary disabled:bg-fill-f2 disabled:hover:opacity-100",
  secondary:
    "bg-fill-f1 text-label-primary hover:bg-fill-f2 disabled:text-label-quaternary disabled:hover:bg-fill-f1",
  outline:
    "bg-transparent border border-sep-s1 text-label-primary hover:bg-fill-f1 disabled:text-label-quaternary disabled:hover:bg-transparent",
  dangerText:
    "bg-fill-f1 text-danger hover:bg-light-red-bg disabled:text-label-quaternary disabled:hover:bg-fill-f1",
};

export function KimiButton({
  variant = "secondary",
  size = "default",
  block = false,
  className = "",
  ...rest
}: KimiButtonProps) {
  const sizeCls =
    size === "hero"
      ? "h-11 rounded-lg px-6 text-base font-medium"
      : "h-8 rounded-md px-3 text-sm font-medium";
  return (
    <button
      type="button"
      className={[
        "inline-flex items-center justify-center gap-1 whitespace-nowrap transition-colors duration-150",
        "active:scale-[0.97] active:transition-transform active:duration-100",
        "disabled:cursor-not-allowed disabled:active:scale-100",
        sizeCls,
        variantCls[variant],
        block ? "w-full" : "",
        className,
      ].join(" ")}
      {...rest}
    />
  );
}

// --------------------------------------------------------------------------
// 卡片与标题（§2/§3）
// --------------------------------------------------------------------------
export function Card({ children, className = "" }: { children: ReactNode; className?: string }) {
  return (
    <section className={`rounded-lg border border-sep-s1 bg-white p-5 ${className}`}>
      {children}
    </section>
  );
}

export function CardTitle({ children, hint }: { children: ReactNode; hint?: string }) {
  return (
    <div className="mb-4">
      <h2 className="text-base font-medium leading-6 text-label-primary">{children}</h2>
      {hint ? <p className="mt-1 text-xs leading-[18px] text-label-tertiary">{hint}</p> : null}
    </div>
  );
}

/** 表单字段标签（14px 次级文字） */
export function FieldLabel({ children, optional }: { children: ReactNode; optional?: boolean }) {
  return (
    <label className="mb-1.5 block text-sm leading-5 text-label-secondary">
      {children}
      {optional ? <span className="ml-1 text-xs text-label-tertiary">（可选）</span> : null}
    </label>
  );
}

// --------------------------------------------------------------------------
// 提示条（§7）：成功/错误/信息/警告
// --------------------------------------------------------------------------
type AlertKind = "success" | "error" | "info" | "warning";

const alertCls: Record<AlertKind, string> = {
  success: "bg-light-green-bg text-positive",
  error: "bg-light-red-bg text-danger",
  info: "bg-light-blue-bg text-label-primary",
  warning: "bg-light-orange-bg text-label-primary",
};

export function Alert({ kind, children }: { kind: AlertKind; children: ReactNode }) {
  return (
    <div className={`rounded-md px-3 py-2 text-sm leading-5 ${alertCls[kind]}`}>{children}</div>
  );
}

// --------------------------------------------------------------------------
// Segmented Control（§6）
// --------------------------------------------------------------------------
export function Segmented<T extends string>({
  options,
  value,
  onChange,
}: {
  options: Array<{ value: T; label: string }>;
  value: T;
  onChange: (v: T) => void;
}) {
  return (
    <div
      role="tablist"
      className="flex h-10 items-center rounded-xl bg-fill-f2 p-1"
    >
      {options.map((opt) => (
        <button
          key={opt.value}
          role="tab"
          aria-selected={value === opt.value}
          type="button"
          onClick={() => onChange(opt.value)}
          className={[
            "h-8 flex-1 rounded-lg px-3 text-sm text-label-primary transition-all duration-150 ease-out",
            value === opt.value ? "bg-white shadow-[0_1px_3px_rgba(0,0,0,0.08)]" : "bg-transparent hover:bg-fill-f1",
          ].join(" ")}
        >
          {opt.label}
        </button>
      ))}
    </div>
  );
}

// --------------------------------------------------------------------------
// 下拉（§7）
// --------------------------------------------------------------------------
export const KimiSelect = forwardRef<HTMLSelectElement, SelectHTMLAttributes<HTMLSelectElement>>(
  function KimiSelect({ className = "", ...rest }, ref) {
    return <select ref={ref} className={`kimi-field ${className}`} {...rest} />;
  },
);

// --------------------------------------------------------------------------
// 数据表（§7）：表头 12px 三级文字 + 下分隔线；行 hover f1；数字右对齐
// --------------------------------------------------------------------------
export interface TableCol {
  key: string;
  title: ReactNode;
  numeric?: boolean;
  render?: (row: Record<string, unknown>, idx: number) => ReactNode;
}

export function DataTable({
  cols,
  rows,
  maxHeight,
}: {
  cols: TableCol[];
  rows: Array<Record<string, unknown>>;
  maxHeight?: number;
}) {
  return (
    <div
      className="overflow-auto rounded-md border border-sep-s1"
      style={maxHeight ? { maxHeight } : undefined}
    >
      <table className="w-full border-collapse text-sm">
        <thead className="sticky top-0 bg-white">
          <tr>
            {cols.map((c) => (
              <th
                key={c.key}
                className={`whitespace-nowrap border-b border-sep-s1 px-3 py-2 text-xs font-normal leading-[18px] text-label-tertiary ${
                  c.numeric ? "text-right" : "text-left"
                }`}
              >
                {c.title}
              </th>
            ))}
          </tr>
        </thead>
        <tbody>
          {rows.map((row, i) => (
            <tr key={i} className="transition-colors duration-150 hover:bg-fill-f1">
              {cols.map((c) => (
                <td
                  key={c.key}
                  className={`whitespace-nowrap px-3 py-1.5 leading-5 text-label-primary ${
                    c.numeric ? "text-right tabular-nums" : "text-left"
                  }`}
                >
                  {c.render ? c.render(row, i) : String(row[c.key] ?? "")}
                </td>
              ))}
            </tr>
          ))}
          {rows.length === 0 ? (
            <tr>
              <td colSpan={cols.length} className="px-3 py-4 text-center text-label-tertiary">
                暂无数据
              </td>
            </tr>
          ) : null}
        </tbody>
      </table>
    </div>
  );
}
