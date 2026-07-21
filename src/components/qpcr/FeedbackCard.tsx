/**
 * FeedbackCard.tsx — 意见反馈
 *
 * 通过 Web3Forms（https://web3forms.com）把表单内容直接发送到站长邮箱：
 * 站长邮箱不会出现在前端代码里，访问者看到的是表单而不是邮箱地址。
 * Access Key 是公开的表单端点标识（非密钥、不含邮箱），可安全提交到仓库。
 *
 * 配置方法（一次性）：
 *   1. 打开 https://web3forms.com ，输入站长邮箱，点 Create Access Key；
 *   2. 到邮箱查收 Access Key；
 *   3. 填到下面常量里，重新构建部署即可。
 */
import { useState, type FormEvent } from "react";
import { Alert, Card, CardTitle, FieldLabel, KimiButton, KimiSelect } from "./ui.tsx";

const WEB3FORMS_ACCESS_KEY = "f909e50c-a29d-414f-ad81-9e629e6fab7c"; // Web3Forms Access Key（公开表单标识，不含邮箱）

const FEEDBACK_TYPES = ["问题反馈", "功能建议", "仪器格式适配请求", "其他"];

type Status = "idle" | "sending" | "success" | "error";

export function FeedbackCard() {
  const [type, setType] = useState(FEEDBACK_TYPES[0]);
  const [contact, setContact] = useState("");
  const [message, setMessage] = useState("");
  const [status, setStatus] = useState<Status>("idle");

  const configured = WEB3FORMS_ACCESS_KEY.length > 0;

  const submit = (e: FormEvent) => {
    e.preventDefault();
    if (!message.trim() || status === "sending") return;
    setStatus("sending");
    fetch("https://api.web3forms.com/submit", {
      method: "POST",
      headers: { "Content-Type": "application/json", Accept: "application/json" },
      body: JSON.stringify({
        access_key: WEB3FORMS_ACCESS_KEY,
        subject: `qPCR 计算器反馈｜${type}`,
        from_name: "qPCR ΔΔCt 计算器",
        botcheck: "", // 反垃圾蜜罐（隐藏字段，正常用户为空）
        反馈类型: type,
        联系方式: contact.trim() || "（未留）",
        message: message.trim(),
      }),
    })
      .then((r) => r.json())
      .then((data: { success?: boolean }) => {
        if (data.success) {
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
