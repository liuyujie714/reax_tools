/* Classic Web Worker that runs the ReaxTools WASM core in the background.
 * The Emscripten modules are UMD-style scripts, so we load them with
 * importScripts(); the threaded build can then spawn its own pthread workers.
 */

let modulePromise = null;

function readDirRecursive(FS, files) {
  const walk = (dir) => {
    for (const name of FS.readdir(dir)) {
      if (name === "." || name === "..") continue;
      const full = dir === "/" ? `/${name}` : `${dir}/${name}`;
      const stat = FS.stat(full);
      if (FS.isDir(stat.mode)) {
        walk(full);
      } else {
        files[full.replace(/^\/out\//, "")] = FS.readFile(full, { encoding: "utf8" });
      }
    }
  };
  walk("/out");
}

self.onmessage = async (event) => {
  const msg = event.data;
  if (!msg || msg.type !== "run") return;
  const logs = [];
  const postLog = (text) => {
    const line = String(text);
    logs.push(line);
    self.postMessage({ type: "log", line });
  };
  self.postMessage({ type: "started" });
  try {
    if (!modulePromise) {
      importScripts(msg.moduleUrl);
      modulePromise = createReaxTools({
        print: postLog,
        printErr: postLog,
        locateFile: (file) => {
          const base = msg.moduleUrl.slice(0, msg.moduleUrl.lastIndexOf("/") + 1);
          return base + file;
        },
        mainScriptUrlOrBlob: msg.moduleUrl,
      });
    }
    const mod = await modulePromise;
    if (!mod.FS.analyzePath("/" + msg.inputName).exists) {
      mod.FS.writeFile(msg.inputName, new Uint8Array(msg.inputBytes));
    }

    const args = msg.args;
    const ptrs = args.map((s) => {
      const p = mod._malloc(s.length + 1);
      mod.stringToUTF8(s, p, s.length + 1);
      return p;
    });
    const argvPtr = mod._malloc(ptrs.length * 4);
    ptrs.forEach((p, i) => mod.setValue(argvPtr + i * 4, p, "i32"));

    const t0 = Date.now();
    const ret = mod.ccall("reax_run_analysis", "number", ["number", "number"], [args.length, argvPtr]);
    const elapsedMs = Date.now() - t0;
    ptrs.forEach((p) => mod._free(p));
    mod._free(argvPtr);

    const files = {};
    readDirRecursive(mod.FS, files);
    self.postMessage({ type: "done", ret, files, logs, elapsedMs });
  } catch (error) {
    self.postMessage({
      type: "done",
      ret: -1,
      files: {},
      logs,
      elapsedMs: 0,
      error: String(error && error.stack ? error.stack : error),
    });
  }
};
