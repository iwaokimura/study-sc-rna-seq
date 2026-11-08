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

## 解説の方針

- 定義・数式をきちんと示す．確率モデルやアルゴリズムは定義と数式で述べ，仮定を明示する．直観的説明は補足にとどめる．
- 生物学的事実は，数理的に必要な範囲で簡潔に導入する．
- 生物学的な解釈には，データや手法の限界（例：PBMC は分化系列データではない）を併記する．

## 表記ルール

- 日本語の文章は **だ・である調**．句読点は **全角「，」「．」** を用いる（「、」「。」は使わない）．
- Markdown 中の数式は LaTeX 記法で書く．インラインは `$...$`，ディスプレイは `$$...$$` で囲む（例：`$\tilde{x}_{ij} = x_{ij} / s_i$`，`$\mathbb{Z}_{\geq 0}^{m \times n}$`）．Unicode の数学記号やプレーンテキストでの数式表記は避ける．
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

- Claude がセルを編集・実行してよい．実行は `course/` 直下の仮想環境で行う：

  ```bash
  cd course
  source .venv/bin/activate      # 未作成なら: uv venv .venv && uv pip install -r requirements.txt
  jupyter nbconvert --to notebook --execute --inplace sc-rna-seq.ipynb
  ```

- データの相対パスは `course/` 基準（`data/filtered_gene_bc_matrices/hg19/`）．
- 出力（図・実行結果）はノートブックに含めたままコミットする．
- パッケージを追加したら `uv pip freeze > requirements.txt` で更新する．

## Git

- コミットはユーザーから依頼されたときのみ，`main` に行う．
- push はしない．
- コミットメッセージに Claude の session ID（セッション URL を含む）を入れない．
