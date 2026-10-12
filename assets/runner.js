// Runs one of the KOALA tools in the browser.
//
// Pyodide (Python compiled to WebAssembly) runs the unchanged tool script;
// its page is shown here, and its form is answered here instead of by a
// local server. Nothing is sent anywhere: all the maths runs on your machine.

const PYODIDE_BASE = "https://cdn.jsdelivr.net/npm/pyodide@314.0.7/";

const TOOLS = [
  ["mosaic.html", "Mosaic planner"],
  ["offsets.html", "Small offsets"],
  ["positions.html", "Offset between positions"],
  ["focus.html", "Focus"],
];

const tool = document.body.dataset.tool;
let runTool = null;

function bar() {
  const here = location.pathname.split("/").pop();
  const links = TOOLS.map(([href, label]) =>
    `<a href="${href}" style="color:${href === here ? "#fff" : "#9fb3d9"};` +
    `text-decoration:none;margin-right:16px;${href === here ?
      "font-weight:600;" : ""}">${label}</a>`).join("");
  return `<div style="font:13px -apple-system,'Segoe UI',Helvetica,Arial,` +
    `sans-serif;background:#0b1020;color:#9fb3d9;padding:8px 24px;` +
    `display:flex;flex-wrap:wrap;align-items:center;gap:4px">` +
    `<a href="./" style="color:#fff;text-decoration:none;font-weight:700;` +
    `margin-right:22px">&#9733; KOALA tools</a>${links}` +
    `<a href="tools/${tool}.py" download style="margin-left:auto;` +
    `color:#9fb3d9">Download ${tool}.py</a></div>`;
}

function show(query) {
  let page = runTool(tool, query);
  page = page.replace(/<form method="get" action="\/">/g,
                      '<form method="get" action="">');
  const doc = new DOMParser().parseFromString(page, "text/html");
  const y = window.scrollY;
  document.title = doc.title;
  document.head.innerHTML =
    '<meta charset="utf-8"><meta name="viewport" ' +
    'content="width=device-width, initial-scale=1">' + doc.head.innerHTML;
  document.body.innerHTML = bar() + doc.body.innerHTML;
  // stretch the bar over the tool page's own margin
  const cs = getComputedStyle(document.body);
  const nav = document.body.firstElementChild;
  nav.style.margin = `-${cs.marginTop} -${cs.marginRight} ` +
    `${cs.marginTop} -${cs.marginLeft}`;
  window.scrollTo(0, y);
  const first = document.querySelector("[autofocus]");
  if (first && !y) first.focus({ preventScroll: true });
}

// The tool's form: run the tool again instead of asking a server.
document.addEventListener("submit", (event) => {
  if (!runTool) return;
  event.preventDefault();
  const data = new FormData(event.target, event.submitter);
  const query = new URLSearchParams(data).toString();
  history.pushState(null, "", "?" + query);
  try {
    show(query);
  } catch (err) {
    failed(err);
  }
});

window.addEventListener("popstate", () => {
  if (runTool) show(location.search.slice(1));
});

function status(text) {
  const el = document.getElementById("boot-status");
  if (el) el.textContent = text;
}

function failed(err) {
  console.error(err);
  document.body.innerHTML = bar() +
    '<div style="font:15px sans-serif;padding:40px">' +
    "<h2>Sorry, the tool could not start in this browser.</h2>" +
    "<p>You can always run it on your own computer: download " +
    `<a href="tools/${tool}.py" download>${tool}.py</a> and run ` +
    `<code>python3 ${tool}.py</code> (Python 3.8 or newer, nothing else ` +
    "to install).</p><pre style='color:#a00'>" +
    String(err).replace(/</g, "&lt;") + "</pre></div>";
}

async function boot() {
  try {
    status("Loading Python into your browser...");
    const py = await loadPyodide({ indexURL: PYODIDE_BASE });
    status("Loading the tool...");
    for (const [file, url] of [[tool + ".py", `tools/${tool}.py`],
                               ["bridge.py", "assets/bridge.py"]]) {
      const reply = await fetch(url);
      if (!reply.ok) throw new Error(`Could not load ${url}`);
      py.FS.writeFile("/home/pyodide/" + file, await reply.text());
    }
    py.runPython("import sys; sys.path.insert(0, '/home/pyodide')");
    runTool = py.runPython("from bridge import run_tool; run_tool");
    show(location.search.slice(1));
  } catch (err) {
    failed(err);
  }
}

boot();
