# CLAUDE.md

## このリポジトリについて

数学者向けの scRNA-seq 解析入門ハンズオンの学習リポジトリ．Scanpy（Python）と 10x Genomics PBMC 3k データで，QC から軌跡推定までの典型的パイプラインを体験する．概要・構成・環境構築手順は [README.md](README.md) を参照．

## ユーザーについて

- 数学者．分子生物学・医学生理学は初学者．
- Python には慣れている（numpy / pandas などの基礎説明は不要）．
- 統計・機械学習にも詳しい（検定・PCA・クラスタリングなどの初歩的説明は不要）．

## Claude の役割

1. **概念の解説**：生物学用語・解析手法を数学の言葉で説明する．
2. **ノートブック実装**：`course/sc-rna-seq.ipynb` のセル追加・修正と実行確認．
3. **教材の執筆・改訂**：`course_plan/` のノートや README の加筆・整理．
4. **文献・事実の検証**：参考文献，パラメータの妥当性，ツールの仕様を一次資料（論文・公式ドキュメント）で確認する．記憶に頼った断定は避け，未確認の事項はその旨を明記する．

## 文献調査

文献の検索には次のツールを使う：

- **数学の文献**：zbMATH Open の MCP（`mcp__zbmath__*`）．
- **より一般的な文献**：CiNii の MCP（`mcp__cinii__*`）．
- **上記以外**（生物学・バイオインフォマティクスの論文など）：Perplexity の MCP（`mcp__perplexity__*`），OpenAlex，一般の Web 検索を使ってよい．

## 解説の方針

- 定義・数式をきちんと示す．確率モデルやアルゴリズムは定義と数式で述べ，仮定を明示する．直観的説明は補足にとどめる．
- 生物学的事実は，数理的に必要な範囲で簡潔に導入する．
- 生物学的な解釈には，データや手法の限界（例：PBMC は分化系列データではない）を併記する．

## 表記ルール

- 日本語の文章は **だ・である調**．句読点は **全角「，」「．」** を用いる（「、」「。」は使わない）．
- Markdown 中の数式は LaTeX 記法で書く．インラインは `$...$`，ディスプレイは `$$...$$` で囲む（例：`$\tilde{x}_{ij} = x_{ij} / s_i$`，`$\mathbb{Z}_{\geq 0}^{m \times n}$`）．Unicode の数学記号やプレーンテキストでの数式表記は避ける．
- Markdown（GitHub で表示するもの）では，インライン数式の `$` が全角の句読点・括弧（「，」「．」「（」「）」など）に隣接する場合，その間に半角スペースを入れる（GitHub で数式として認識させるため）．LaTeX 文書の規則とは異なる．
  - 例：`行列（ $m$ = 遺伝子数），` `ライブラリサイズ， $s_i = \sum_j x_{ij}$ ）は`
- `course_plan/scrna_hands_on_course.md` の参考文献は `\[Ref. n\]` 形式で引用し，末尾の参考文献リストに DOI / URL 付きで追加する．

## LaTeX 文書

日本語の LaTeX 文書を作成する場合：

- `lualatex` でコンパイルできるものにする．
- 日本語化は `luatexja-preset` を `haranoaji`，`deluxe` オプション付きで用いる：

  ```latex
  \usepackage[haranoaji,deluxe]{luatexja-preset}
  ```

- 句読点は全角「，」「．」．
- インライン数式のデリミタ `$` の前後にスペースを入れない．
  - 良い例：`リーマンゼータ関数$\zeta(s)$を次のように定義する．`
  - 悪い例：`リーマンゼータ関数 $\zeta(s)$ を次のように定義する．`

## ファイルの扱い

- `course_plan/scrna_hands_on_course.md`（詳細版）が正．内容の更新はこちらのみに行う．
- `course_plan/sc-RNAseq-analysis.md`（初版）は記録として残し，編集しない．
- `course/data/` は git 管理外．データは README の手順で取得する．

## ノートブック

- Claude がセルを編集・実行してよい．実行はリポジトリ直下の uv 環境で行う：

  ```bash
  uv sync    # 環境が未作成・古い場合
  uv run jupyter nbconvert --to notebook --execute --inplace --allow-errors --ExecutePreprocessor.timeout=1800 course/sc-rna-seq.ipynb
  ```

- Python は 3.13（リポジトリ直下の `.python-version`）．
- 全体の実行には約 10 分かかる（大半はモジュール 7 の RNA velocity）．
- モジュール 6 発展版（`course/module6-multidonor.ipynb`，約 1 分）は依存グループ `multidonor` が必要で，既定の `uv sync` を実行すると PyDESeq2 が外れる．次のように実行する：

  ```bash
  uv run --group multidonor jupyter nbconvert --to notebook --execute --inplace --allow-errors --ExecutePreprocessor.timeout=1800 course/module6-multidonor.ipynb
  ```

- データの相対パスは `course/` 基準（`data/filtered_gene_bc_matrices/hg19/`）．
- 出力（図・実行結果）はノートブックに含めたままコミットする．
- パッケージは `uv add <pkg>`（発展版用など任意のものは `uv add --group <名前> <pkg>`）で追加し，`pyproject.toml` と `uv.lock` をコミットする．`uv pip install` で直接入れない．

## Git

- コミットはユーザーから依頼されたときのみ，`main` に行う．
- push はしない．
- コミットメッセージに Claude の session ID（セッション URL を含む）を入れない．
