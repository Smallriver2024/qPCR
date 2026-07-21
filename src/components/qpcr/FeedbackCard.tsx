/**
 * FeedbackCard.tsx — 意见反馈
 *
 * 通过 FormSubmit（https://formsubmit.co）把表单内容直接发送到站长邮箱：
 * 前端代码里只放「随机字符串别名」，站长邮箱不出现在代码和页面上，
 * 访问者看到的是表单而不是邮箱地址。AJAX 端点支持跨域，静态站可用。
 *
 * 配置方法（一次性）：
 *   1. 向 https://formsubmit.co/ajax/{站长邮箱} POST 一次，触发激活邮件；
 *   2. 到邮箱点击 "Activate Form" 完成激活；
 *   3. 确认邮箱地址后 FormSubmit 会提供随机字符串（Invisible emails 功能），
 *      填到下面常量里，重新构建部署即可。
 */
import { useState, type FormEvent } from "react";
import { Alert, Card, CardTitle, FieldLabel, KimiButton, KimiSelect } from "./ui.tsx";

const FORMSUBMIT_ALIAS = ""; // ← 在这里填入 FormSubmit 随机字符串别名（不含邮箱）

const FEEDBACK_TYPES = ["问题反馈", "功能建议", "仪器格式适配请求", "其他"];

type Status = "idle" | "sending" | "success" | "error";

export function FeedbackCard() {
  const [type, setType] = useState(FEEDBACK_TYPES[0]);
  const [contact, setContact] = useState("");
  const [message, setMessage] = useState("");
  const [status, setStatus] = useState<Status>("idle");

  const configured = FORMSUBMIT_ALIAS.length > 0;

  const submit = (e: FormEvent) => {
    e.preventDefault();
    if (!message.trim() || status === "sending") return;
    setStatus("sending");
    fetch(`https://formsubmit.co/ajax/${FORMSUBMIT_ALIAS}`, {
      method: "POST",
      headers: { "Content-Type": "application/json", Accept: "application/json" },
      body: JSON.stringify({
        _subject: `qPCR 计算器反馈｜${type}`,
        _template: "table", // 邮件用表格模板，更易读
        _captcha: "false", // 站内表单 + 蜜罐即可，不要求访客过验证码
        _honey: "", // 反垃圾蜜罐（隐藏字段，正常用户为空）
        反馈类型: type,
        联系方式: contact.trim() || "（未留）",
        message: message.trim(),
      }),
    })
      .then((r) => r.json())
      .then((data: { success?: boolean | string }) => {
        if (data.success === true || data.success === "true") {
          setStatus("success");
          setMessage("");
        } else {
          setStatus("error");
        }
      })
      .catch(() => setStatus("error"));
  };

  return (
    <Card>
      <CardTitle hint="提交后会直接发送到站长邮箱，无需登录；留下联系方式便于回复">
        意见反馈
      </CardTitle>
      {configured ? (
        <form className="space-y-4" onSubmit={submit}>
          <div className="grid grid-cols-1 gap-4 sm:grid-cols-2">
            <div>
              <FieldLabel>反馈类型</FieldLabel>
              <KimiSelect value={type} onChange={(e) => setType(e.target.value)}>
                {FEEDBACK_TYPES.map((t) => (
                  <option key={t} value={t}>{t}</option>
                ))}
              </KimiSelect>
            </div>
            <div>
              <FieldLabel optional>联系方式（邮箱 / 微信，便于回复）</FieldLabel>
              <input
                value={contact}
                onChange={(e) => setContact(e.target.value)}
                className="kimi-field"
                placeholder="可不填"
              />
            </div>
          </div>
          <div>
            <FieldLabel>反馈内容</FieldLabel>
            <textarea
              value={message}
              onChange={(e) => setMessage(e.target.value)}
              rows={4}
              required
              className="kimi-field w-full text-[13px] leading-5"
              placeholder="描述遇到的问题（最好附上仪器型号与导出方式），或希望增加的功能……"
            />
          </div>
          <div className="flex items-center gap-3">
            <KimiButton
              variant="primary"
              type="submit"
              disabled={!message.trim() || status === "sending"}
            >
              {status === "sending" ? "发送中…" : "提交反馈"}
            </KimiButton>
            {status === "success" ? (
              <span className="text-sm leading-5 text-positive">已发送，感谢反馈！</span>
            ) : null}
          </div>
          {status === "error" ? (
            <Alert kind="error">发送失败，请检查网络后重试；若多次失败，可到 GitHub 仓库提 Issue。</Alert>
          ) : null}
        </form>
      ) : (
        <Alert kind="info">反馈通道配置中，即将开放。</Alert>
      )}
    </Card>
  );
}
