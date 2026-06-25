# ReaxTools v2.1 更新说明

ReaxTools v2.1 是一次比较大的结构重构版本，核心目标是让程序输出更稳定、更干净、更可审计。

## 主要变化

1. 统一入口

现在用户只需要使用一个命令：

```bash
reax_tools -f traj.xyz -o output_dir
reax_tools plot -f output_dir
reax_tools network -f output_dir
reax_tools events -f output_dir
reax_tools counts -f output_dir
```

直接 `reax_tools -f ...` 等价于分析模式；`plot/network/events/counts/flow/focus` 等子命令负责后处理和画图。

2. C++ 和 Python 职责分离

C++ 只负责高性能分析和可审计 raw 输出：

- `species_count.csv`
- `bond_count.csv`
- `atom_bonded_num_count.csv`
- `ring_count.csv`
- `reaction_events.csv`
- `reaction_event_pairs.csv`
- `transfer_flow.csv`
- `molecules.json`
- `reax_tools.log`
- `reax_tools_manifest.json`

Python 负责过滤、画图和展示。这样以后大多数用户体验和绘图改动都可以在 Python 层快速更新，不需要频繁改 C++ 核心。

3. 移除不必要外部依赖

本地默认构建不再依赖 RDKit、Boost、Graphviz 或打包的动态库。安装和分发都更轻。

4. 反应网络审计增强

v2.1 新增 `reaction_event_pairs.csv`，它记录每个 raw reaction event 内 reactant-product 的原子重叠 pair。

`transfer_flow.csv` 是 `reaction_event_pairs.csv` 的逐边聚合结果。也就是说，完整 network 不再只是“画出来的图”，而是可以从 raw reaction events 审计出来的网络。

我们增加了审计工具：

```bash
python3 tools/check_reax_outputs.py output_dir
python3 tools/audit_event_network_consistency.py output_dir
```

在 `energetic_v3` 和 `gpumd_v3` 测试体系中，event-network 一致性审计均通过。

5. network 和 flow 的定位更清楚

- `network`：完整的、可审计的物质转移网络视图，是反应网络分析的核心。
- `flow`：实验性 Sankey 风格叙事图，用来帮助快速看主线，但会压缩和丢弃部分循环/同层信息。

因此，正式分析和可审计结论建议以 raw CSV 和 `network` 为准；`flow` 更适合作为探索和展示辅助。

6. README 和输出风格更新

新版 README 以 `test/energetic_v3` 为样例，展示 species、reaction events、transfer network 和 experimental flow。命令行输出也比旧版更克制，详细信息放到 log 和 manifest 中。

## 推荐安装

```bash
git clone https://github.com/tgraphite/reax_tools
cd reax_tools
bash install_reax_tools.sh
```

## 推荐使用

```bash
reax_tools -f trajectory.xyz -o reax_tools_output
reax_tools plot -f reax_tools_output
reax_tools network -f reax_tools_output
```

LAMMPS dump 文件通常需要指定元素类型：

```bash
reax_tools -f dump.lammpstrj -t C,H,O,N -o reax_tools_output
```

## 一句话总结

ReaxTools v2.1 的重点不是增加花哨功能，而是把反应网络相关结果从“能画图”推进到“可审计、可复现、可长期维护”。完整 network 是核心可信结果；flow 目前仍是实验性展示功能。
