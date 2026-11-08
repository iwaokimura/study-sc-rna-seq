# 計画：モジュール 6 発展版 — 複数ドナーデータによる擬似バルク差次発現解析

> **状態**：計画（未着手）．記録日 2026-11-08（同日改訂）．
> **位置づけ**：モジュール 6 の擬似バルク節（PBMC 3k の無作為分割による「集約の実習」）の続き．初版の計画（同日記録）を，データと API の確認結果とユーザーの決定を反映して改訂した．

## Context

現行モジュール 6 の擬似バルク節は，単一ドナーの PBMC 3k を無作為に分けた仮サンプルで集約 $Y = AX$ を計算するだけで，検定はしていない．発展版では，生物学的反復（ドナー）と条件を持つ実データで，

1. ドナー × 条件 × 細胞型ごとの擬似バルク行列を作り，
2. 負の二項 GLM（PyDESeq2）で条件間の差次発現を検定し，
3. 細胞単位の Wilcoxon 検定（擬似反復）と比べて，偽陽性の増え方を実験で確かめる．

前提の整備（Python 3.13，`pyproject.toml` + `uv.lock`）は完了済み．

**決定事項（ユーザー回答）**
- 既存ノートブックには追記せず，**別ノートブック** `course/module6-multidonor.ipynb` にする．
- pertpy は使わない．データは **curl で直接取得** し，**PyDESeq2 を直接呼ぶ**（NB-GLM の各段階を見える形にする）．

## 確認済みの事実（2026-11-08）

- **データ**：`https://exampledata.scverse.org/pertpy/kang_2018.h5ad`（38,356,412 バイト，HTTP 200）．pertpy の `kang_2018()` もこの URL から取得している（pertpy のソースで確認）．sc-best-practices の章で使われているデータと形が同じ．
  - 24,673 細胞 × 15,706 遺伝子．`X` は整数値の生カウント（CSC 疎行列，最大値 3828）．`layers` は空．
  - `obs`：`label`（`ctrl` 12,315 / `stim` 12,358），`replicate`（`patient_101` など 8 名），`cell_type`（8 種：CD4 T cells, CD14+ Monocytes, B cells, NK cells, CD8 T cells, FCGR3A+ Monocytes, Dendritic cells, Megakaryocytes）．
  - CD14+ Monocytes はドナー × 条件の 16 サンプルすべてで 126〜859 細胞ある．Megakaryocytes は 0 細胞のサンプルがあり，CD8 T cells は最少 8 細胞．
- **sc-best-practices の章**（https://www.sc-best-practices.org/conditions/differential-gene-expression/ ）：
  - 前処理は `filter_cells(min_genes=200)` と `filter_genes(min_cells=3)`．
  - 擬似バルクは decoupler の `pseudobulk(sample_col="sample", groups_col="cell_type", mode="sum")`．
  - サンプルの除外条件は `min_cells=10, min_counts=1000`．これは「よく使われる目安」とされている．
  - 遺伝子の除外は `filter_by_expr`（`min_count=10, min_total_count=15` など）．
  - 主解析は CD14+ Monocytes で，計画式は `~ label`．使用版は pydeseq2 0.5.4．
  - 擬似反復の文献として Squair et al. 2021（https://doi.org/10.1038/s41467-021-25960-2 ），Zimmerman et al. 2021（https://doi.org/10.1038/s41467-021-21038-1 ），Junttila et al. 2022 を挙げる．
- **PyDESeq2 0.5.x の API**：
  - `DeseqDataSet(counts=..., metadata=..., design="~replicate + label", refit_cooks=True, inference=DefaultInference(n_cpus=...))` を作り，`dds.deseq2()` を呼ぶ．
  - `DeseqStats(dds, contrast=["label", "stim", "ctrl"], inference=...)` を作り，`.summary()` を呼ぶと結果が `.results_df` に入る．
  - LFC の縮小推定は `.lfc_shrink(coeff=...)`．
  - `design_factors` は非推奨で，計画式（formula）が推奨されている．
  - 正確な引数名と `results_df` の列は，実装時に導入した版で確認する（Perplexity の回答で，公式ドキュメント v0.5.1〜0.5.4 を出典とする．一次資料での最終確認は実装時）．

## 実装方針

### 1. パッケージ

```bash
uv add --group multidonor pydeseq2
```

- `pyproject.toml` の `[dependency-groups]` と `uv.lock` を更新してコミットする．
- uv の既定では，`uv sync` で全グループが入る（入れたくない場合は `--no-group multidonor`）．この挙動は実装時に確認し，README に書く．
- 既存の主パッケージの版は動かさない．ロック後に `git diff uv.lock` で確認する．

### 2. データ取得（README に追記）

```bash
cd course
curl -L https://exampledata.scverse.org/pertpy/kang_2018.h5ad -o data/kang_2018.h5ad
```

`course/data/` は git 管理外のまま．

### 3. ノートブック `course/module6-multidonor.ipynb` の構成

既存ノートブックの書き方（日本語 markdown の見出しと注釈，`import scanpy as sc`，データパスは `course/` 基準）に合わせる．

1. **導入**（markdown）：PBMC 3k との違い，データの由来（ループス患者 8 名の PBMC．各ドナーを未処置と IFN-β 刺激に分けた対応のある設計）．
2. **読み込みと確認**：
   - `sc.read_h5ad('data/kang_2018.h5ad')` で読み込む．
   - 生カウントであることを `assert` で確かめる．
   - `pd.crosstab` でドナー × 条件 × 細胞型の細胞数を表示する．
   - sc-best-practices と同じ最小限の QC（`min_genes=200`，`min_cells=3`）を行う．
3. **擬似バルク化**：
   - 所属行列 $A \in \{0,1\}^{G \times n}$ による定式化 $Y = AX$ を，本編モジュール 6 から再掲する．
   - 実装は `sc.get.aggregate(adata, by=['replicate','label','cell_type'], func='sum')`．本編で $A$ との一致を確かめた関数を再利用する．
   - 細胞数 10 未満，または総カウント 1000 未満のサンプルを除外する（sc-best-practices に合わせる）．除外されたサンプルを表示する．
4. **探索**（短く）：擬似バルクを CPM と $\log(1+x)$ で変換して PCA にかけ，条件・ドナー・細胞型で色分けする．
5. **NB-GLM の定式化**（markdown，主眼）：CD14+ Monocytes に絞り，遺伝子 $j$ ，擬似バルクサンプル $g$ について
   $$Y_{gj} \sim \mathrm{NB}(\mu_{gj}, \alpha_j), \qquad \mathrm{Var}(Y_{gj}) = \mu_{gj} + \alpha_j \mu_{gj}^2, \qquad \log \mu_{gj} = \log s_g + \mathbf{x}_g^\top \boldsymbol{\beta}_j$$
   とおく．次を数式で説明する：
   - サイズ因子 $s_g$ の median-of-ratios 推定．
   - 分散 $\alpha_j$ の遺伝子ごとの推定，平均–分散トレンドへの経験ベイズ的な縮小（MAP 推定）．
   - 係数の Wald 検定．
   - Benjamini–Hochberg による多重検定補正．
6. **遺伝子フィルタ**：CD14+ Monocytes のサンプルで，十分に発現する遺伝子だけを残す．decoupler は使わず，numpy で簡単な規則を書く．たとえば「CPM がある閾値以上のサンプルが，少ない方の群のサンプル数以上ある」という edgeR の filterByExpr に倣った規則にし，規則は数式で明示する．
7. **PyDESeq2 による検定**：対応のある設計 `~ replicate + label` を主にする（ドナー効果を固定効果として除く）．比較として `~ label`（sc-best-practices の計画式）も実行し，次を表にする：
   - 有意遺伝子数（ $\mathrm{padj} < 0.05$ ）．
   - 上位遺伝子．ISG15 などインターフェロン応答遺伝子が上位に来るはずなので，それを確認する．
   - volcano plot．
8. **擬似反復の実験**（このノートブックの山場）：
   - (a) **実データの比較**：同じ CD14+ Monocytes の細胞で，`sc.tl.rank_genes_groups(groupby='label', method='wilcoxon')` を実行し，有意遺伝子数を擬似バルクの結果と比べる．
   - (b) **帰無データでの偽陽性**：`ctrl` の細胞だけを使い，8 ドナーを無作為に 4 対 4 の仮の群に分ける（真の差はない）．次の 2 つで有意遺伝子数を数える：
     - 細胞単位の Wilcoxon．
     - 擬似バルク PyDESeq2（`~ group`）．
     数回（例：分け方 5 通り，乱数シードは固定）繰り返し，前者で偽陽性が多いことを表で示す．これは Squair et al. 2021 の主張の縮小版である．計算時間を見て回数を決める．
9. **まとめと限界**（markdown）：
   - 固定効果モデルとしての扱いと，混合モデルとの関係．
   - ドナー数 8 という少なさ．
   - 細胞型ラベルは配布元の注釈をそのまま使っていること．

### 4. 文書の更新

- **`course_plan/scrna_hands_on_course.md`**：モジュール 6 に「発展：複数ドナーでの擬似バルク DE」の節を追加する．中身は上記 3–5，7–8 の数式と要点．ノートブックへの参照を付ける．
- **同じファイルの参考文献**：次を `\[Ref. n\]` 形式で追加し，書誌は Crossref などで確認する．
  - Kang et al. 2018（https://doi.org/10.1038/nbt.4042 ）
  - Love et al. 2014 DESeq2（https://doi.org/10.1186/s13059-014-0550-8 ）
  - PyDESeq2 の論文（Muzellec et al. 2023, *Bioinformatics*．DOI は実装時に確認）
  - Squair et al. 2021
  - Zimmerman et al. 2021
- **`README.md`**：
  - 構成図に新しいノートブックを追加する．
  - データ取得に Kang データを追加する．
  - 実行例に `uv run jupyter lab course/module6-multidonor.ipynb` を追加する．
  - 注意点の擬似バルクの項に，発展版への案内を付ける．
- **`CLAUDE.md`**：ノートブックの節に，2 冊目のノートブックの実行コマンドと実行時間の目安を追記する．
- **`course_plan/plan_module6_multidonor.md`**：状態を「実施済み」にする．

## 検証

1. `uv lock --check` が通ること．`git diff uv.lock` で既存パッケージの版が変わっていないこと．
2. 次のコマンドでエラーなく完走すること．
   ```bash
   uv run jupyter nbconvert --to notebook --execute --inplace --allow-errors --ExecutePreprocessor.timeout=1800 course/module6-multidonor.ipynb
   ```
   出力にエラーや想定外の警告がないこと．
3. 数値の健全性を確かめる．
   - 擬似バルクの総和が元のカウント総和と一致すること（除外前，`assert`）．
   - CD14+ Monocytes の 16 サンプルがすべて残ること．
   - `stim` 対 `ctrl` の上位にインターフェロン応答遺伝子（ISG15 など）が来ること．
   - 帰無実験で，擬似バルクの偽陽性数が 0 に近く，Wilcoxon の偽陽性数がそれより明らかに多いこと．
4. 既存ノートブック `course/sc-rna-seq.ipynb` が影響を受けないこと（依存グループの追加のみ．必要なら再実行して出力を比べる）．
5. 実行時間を計測し，CLAUDE.md と README に記す．

## コミット

- 計画書：`course_plan/plan_module6_multidonor.md` を上書きして 1 コミット（今回）．
- 実装：依頼時に，パッケージ追加・ノートブック・文書をまとめて，または段階ごとにコミットする．
