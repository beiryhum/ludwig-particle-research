# AGENTS.md

## 1. Repository Mission

本仓库是基于官方 Ludwig 修改的科研求解器。

主要用途：

- 颗粒解析 LBM；
- 颗粒运动与流固耦合；
- 颗粒总力和力矩；
- 压力与黏性贡献；
- 颗粒表面牵引力；
- 二维和三维验证；
- MPI 并行计算。

本仓库只负责求解器代码、测试和求解器级验证。

论文规划、Run Registry 和正式论文数据由科研管理仓库负责。

## 2. Repository Remotes

研究者仓库：

origin
https://github.com/beiryhum/ludwig-particle-research.git

官方上游：

upstream
https://github.com/ludwig-cf/ludwig.git

不得向 upstream 推送研究者修改。

## 3. Required Documents

进行工作前应优先阅读：

- UPSTREAM_VERSION.md
- BUILD_ENVIRONMENT.md
- RESEARCH_MODIFICATIONS.md
- BASELINE_VALIDATION.md
- README.md

这些文件是本仓库的事实来源。

不得用推测替代其中缺失的信息。

## 4. Current Priority

当前工作顺序：

1. 记录真实构建和运行环境；
2. 审计程序入口和时间步流程；
3. 审计颗粒初始化；
4. 审计流固边界动量交换；
5. 审计颗粒力和力矩；
6. 审计 MPI reduction；
7. 审计位置、速度和角速度更新；
8. 固化二维回归；
9. 验证 MPI 一致性；
10. 验证三维固定球；
11. 验证表面牵引力。

在这些步骤完成前，不进行湍流颗粒生产计算。

## 5. Hard Rules

未经明确批准，不得：

- 修改核心求解器代码；
- 大规模重构；
- 改变物理模型；
- 改变坐标或符号；
- 更新回归参考值以掩盖失败；
- 修改参数以暂时隐藏发散；
- 删除原始测试；
- 覆盖已有基准结果；
- 提交高成本 HPC 作业；
- 将大型输出提交到 Git；
- 将密码、Token 或集群凭证写入仓库。

## 6. Unknown Information

不得猜测：

- 表面法向方向；
- 力是流体作用于颗粒还是相反；
- 力矩参考点；
- 压力和黏性贡献的符号；
- MPI 汇总是否重复；
- 当前输出是否为局部量或全局量；
- 研究者修改过的完整文件列表。

无法通过代码和文档确认时，写：

- [UNKNOWN]
- [TO VERIFY]

## 7. Read-Only Audit Requirements

只读审计必须输出：

1. 程序入口；
2. 时间步调用流程；
3. 颗粒初始化；
4. 流固边界处理；
5. 颗粒力累加；
6. 颗粒力矩累加；
7. MPI reduction；
8. 位置更新；
9. 平动速度更新；
10. 角速度更新；
11. 压力、速度和应力输出；
12. 表面牵引力状态。

每个判断必须给出：

- 文件；
- 函数；
- 变量或数据字段；
- 调用关系；
- 证据；
- 不确定性；
- 风险等级。

## 8. Code Change Protocol

修改前必须先输出：

- Scientific purpose；
- Scope；
- Out of scope；
- Files to modify；
- Functions to modify；
- Data structure changes；
- MPI implications；
- Test plan；
- Acceptance criteria；
- Rollback plan；
- P0/P1/P2 risks。

得到人工批准前不得修改。

## 9. Validation Requirements

任何功能修改至少检查：

- 编译；
- 单元测试；
- 原有回归测试；
- 新功能测试；
- 串行或少核结果；
- 多核结果；
- 输出一致性；
- NaN 和异常值；
- 是否改变原动力学。

表面牵引力功能还必须检查：

- 表面坐标；
- 表面法向；
- 面积权重；
- 压力牵引力；
- 黏性牵引力；
- 合力；
- 力矩；
- MPI 汇总；
- 与独立总力比较；
- 网格收敛。

## 10. Risk Classification

P0：

- 力符号错误；
- 法向错误；
- 面积错误；
- MPI 重复或漏加；
- 力矩参考点错误；
- 输出代码改变动力学；
- 压力和黏性分解错误。

P1：

- 网格不足；
- 插值误差；
- 时间步误差；
- 统计长度不足；
- 边界影响。

P2：

- 路径硬编码；
- 重复代码；
- 文档缺失；
- 输出效率低；
- 文件命名混乱。

发现 P0 风险时停止功能开发或生产运行，先完成验证。

## 11. Build Artifacts

不得提交：

- *.o
- *.a
- *.exe
- a.out
- 本地 config.mk
- core.*
- slurm-*.out
- checkpoint
- 大型二进制结果

不得全局忽略 `*.log`，因为部分回归参考结果可能使用日志文件。

## 12. Modification Report

代码修改完成后必须报告：

- 修改文件；
- 修改函数；
- diff 摘要；
- 构建结果；
- 测试结果；
- 失败测试；
- 输出变化；
- 性能变化；
- 物理风险；
- 已知限制；
- 对应 commit。

不得只报告“修改成功”或“测试通过”。

## 13. Validation Case Placement

求解器级验证放在：

research\validation_cases\

包括：

- settling_2d
- couette_2d
- poiseuille_2d
- fixed_sphere_3d
- shear_sphere_3d

目录中只保存：

- input；
- README；
- 运行脚本；
- 关键参考指标；
- 小型后处理程序。

大型结果不得保存在仓库中。