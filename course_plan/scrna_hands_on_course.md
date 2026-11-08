# scRNA-seq ハンズオンコース：数学者のための入門ガイド

> **対象者**：数理的素養と計算機操作の経験を持つが，分子生物学・医学生理学の背景知識をほとんど持たない数学者  
> **ゴール**：定番ツール（Scanpy / Python）と公開データ（10x Genomics PBMC 3k）を使って scRNA-seq 解析の典型的パイプライン全体を体験する  
> **所要時間**：集中実施で 3〜5 日間（コード実行のみなら 1〜1.5 日，概念理解を含めると 3〜4 日が現実的）

---

## 目次

1. [モジュール 0：準備](#module0)
2. [モジュール 1：カウント行列と品質管理（QC）](#module1)
3. [モジュール 2：正規化と高変動遺伝子選択](#module2)
4. [モジュール 3：次元削減 — PCA と UMAP](#module3)
5. [モジュール 4：クラスタリングとグラフ分割](#module4)
6. [モジュール 5：細胞型アノテーション](#module5)
7. [モジュール 6：差次発現解析（DEG）](#module6)
8. [モジュール 7：発展 — 軌跡推定と擬似時間](#module7)
9. [参考文献](#references)

---

## コース全体構成

| モジュール | タイトル | 所要時間目安 |
|---|---|---|
| 0 | 準備：生物学の最小知識と環境構築 | 0.5 日 |
| 1 | データの正体を知る：カウント行列と QC | 0.5 日 |
| 2 | 正規化・高変動遺伝子選択：統計モデルの話 | 0.5 日 |
| 3 | 次元削減：PCA・UMAP の数理 | 0.5 日 |
| 4 | クラスタリングとグラフ分割 | 0.5 日 |
| 5 | 細胞型アノテーション：マーカー遺伝子と検定 | 0.5〜1 日 |
| 6 | 差次発現解析：負の二項分布モデル | 0.5 日 |
| 7 | 発展：軌跡推定と擬似時間 | 0.5〜1 日 |

---

<a id="module0"></a>
## モジュール 0：準備

### 最小限の生物学知識

数学者が scRNA-seq を理解するために必要な生物学的事実は以下の 3 点に絞られる：

1. **細胞** は生命の基本単位．ヒトの体には約 37 兆個の細胞があり，それぞれが同一のゲノム DNA を持つ．
2. **遺伝子発現**（転写）とは，DNA の情報が mRNA に写し取られるプロセス．mRNA の量が「その遺伝子がどれだけ使われているか」を示す．scRNA-seq はこの mRNA 量を各細胞ごとに計測する技術である．
3. **PBMC**（末梢血単核球）とは血液中の免疫細胞の集合体で，T 細胞（CD4+: 25〜60%，CD8+: 5〜30%），B 細胞（5〜10%，CD19+），NK 細胞（10〜30%，CD56+CD3-），単球（CD14+）などを含む．本コースで用いる PBMC データは，これらの異なる細胞種が混在する状況での **分類問題** として捉えると数学的に理解しやすい．

> **数理的直観**：scRNA-seq の出力は $m \times n$ 行列（ $m$ = 遺伝子数 ≈ 20,000， $n$ = 細胞数 ≈ 3,000）であり，各エントリは非負整数（UMI カウント）である．目標はこの高次元空間の点群を意味のある部分集合（細胞型）に分類することである．

### 使用ツールと公開データ

- **メインツール**：Scanpy（Python）— AnnData オブジェクトをベースにした scRNA-seq 解析パッケージ \[Ref. 1\]
- **データ**：10x Genomics の「PBMC 3k」（健常ドナー由来，2,700 細胞）\[Ref. 2\]
- **参考書**：Single-cell best practices（Theis ラボ, Helmholtz Munich）— 無料オンライン教科書 \[Ref. 3\]
- **ベストプラクティス論文**：Luecken & Theis (2019) \[Ref. 4\]，Heumos et al. (2023) \[Ref. 5\]

### 環境構築

`uv` のプロジェクト管理を使う．依存パッケージはリポジトリ直下の `pyproject.toml` に，全パッケージの正確な版は `uv.lock` に記録されている．Python の版は `.python-version`（3.13）で指定する：

```bash
# uv 本体のインストール（未インストールの場合）
curl -LsSf https://astral.sh/uv/install.sh | sh

# リポジトリ直下で，uv.lock どおりの環境を .venv/ に作る
uv sync

# JupyterLab の起動
uv run jupyter lab course/sc-rna-seq.ipynb
```

JupyterLab での作業を推奨する（Google Colab を使う場合は `uv` 不要で，必要なパッケージを `pip install` で直接入れる）．

### データ取得

```bash
mkdir -p data write
curl https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz \
  -o data/pbmc3k_filtered_gene_bc_matrices.tar.gz
tar -xzf data/pbmc3k_filtered_gene_bc_matrices.tar.gz -C data/
```

あるいは Scanpy の組み込み関数で直接ダウンロードすることもできる：

```python
import scanpy as sc
adata = sc.datasets.pbmc3k()
```

---

<a id="module1"></a>
## モジュール 1：データの正体を知る — カウント行列と品質管理（QC）

### AnnData オブジェクト

Scanpy は AnnData という構造体でデータを管理する \[Ref. 1\]：

```python
import scanpy as sc

adata = sc.read_10x_mtx(
    'data/filtered_gene_bc_matrices/hg19/',
    var_names='gene_symbols'
)
adata.var_names_make_unique()
print(adata)  # AnnData object: 2700 obs × 32738 vars
```

`adata.X` は疎行列（scipy sparse matrix）であり， $2{,}700 \times 32{,}738$ の UMI カウント行列である．エントリの大部分はゼロ（スパース性 > 90%）で，これは **ドロップアウト** と呼ばれる技術的アーティファクトに起因する \[Ref. 4\]．

### QC 指標の計算

```python
# ミトコンドリア遺伝子の識別
adata.var['mt'] = adata.var_names.str.startswith('MT-')

# QC 指標の計算
sc.pp.calculate_qc_metrics(
    adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True
)
```

各細胞に対して計算する主要な QC 指標 \[Ref. 4\]\[Ref. 5\]：

| 指標 | 数理的意味 | フィルタ閾値の目安 |
|---|---|---|
| `n_genes_by_counts` | 非ゼロエントリ数（ $\ell_0$ ノルム的） | 200〜2,500 |
| `total_counts` | 各行の $\ell_1$ ノルム（ライブラリサイズ） | 外れ値除去 |
| `pct_counts_mt` | ミトコンドリア遺伝子カウントの割合 | < 5% |

ミトコンドリア遺伝子の高発現は細胞の「死」の指標，遺伝子数が異常に多い細胞はダブレット（2 細胞が 1 つのバーコードに混入）の疑いがある \[Ref. 5\]．

```python
# 可視化
sc.pl.violin(
    adata,
    ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
    jitter=0.4, multi_panel=True
)
sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt')
sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts')

# フィルタリング
adata = adata[adata.obs.n_genes_by_counts < 2500, :]
adata = adata[adata.obs.pct_counts_mt < 5, :].copy()  # ビューを実体にする
sc.pp.filter_genes(adata, min_cells=3)
```

> **発展**：アンビエント RNA（空の液滴に含まれる細胞外 mRNA の混入）除去には **SoupX** が有効である \[Ref. 6\]．ダブレット検出には **scDblFinder** が推奨される \[Ref. 7\]．本コースでは基本フィルタのみ実施する．

---

<a id="module2"></a>
## モジュール 2：正規化と高変動遺伝子選択

### 正規化：なぜ必要か

各細胞間で捕捉された総 RNA 量（ライブラリサイズ， $s_i = \sum_j x_{ij}$ ）は技術的ばらつきにより異なる．単純なサイズファクター正規化（CPM 正規化相当）を適用する \[Ref. 4\]\[Ref. 5\]：

$$\tilde{x}_{ij} = \frac{x_{ij}}{s_i} \cdot 10{,}000$$

その後対数変換 $\log(1 + \tilde{x}_{ij})$ を施す（`log1p` 変換）\[Ref. 4\]：

```python
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata  # 正規化済みデータを保存
```

より高度な正規化として **SCTransform**（Hafemeister & Satija, 2019）がある \[Ref. 8\]．これは細胞ごとのシーケンス深度を共変量とした一般化線形モデル（GLM）を各遺伝子に対して当てはめ，ピアソン残差を正規化値として使用する．また **scran** 法（Lun et al., 2016）はプーリング戦略により細胞間スケーリングをより頑健に推定する \[Ref. 9\]．

### 高変動遺伝子（HVG）の選択

全 32,738 遺伝子のうち，細胞間で統計的に有意な発現変動を示す遺伝子（高変動遺伝子, Highly Variable Genes: HVG）を選択して以降の解析に用いる \[Ref. 4\]\[Ref. 5\]：

```python
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
print(f"HVG 数: {adata.var.highly_variable.sum()}")  # 通常 2,000〜4,000
adata = adata[:, adata.var.highly_variable].copy()
```

スケーリング（平均 0，分散 1）は PCA の前処理として標準的だが，scRNA-seq ではスパース性を損なうというトレードオフがある \[Ref. 3\]：

```python
# 中心化すると疎でなくなるので，先に密行列にする（2,638 細胞 × HVG なら十分小さい）
adata.X = adata.X.toarray()
sc.pp.scale(adata, max_value=10)
```

---

<a id="module3"></a>
## モジュール 3：次元削減 — PCA と UMAP の数理

### PCA：線形次元削減

正規化済みデータ行列 $X \in \mathbb{R}^{n \times m}$（ $n$ = 細胞， $m$ = HVG 数）に対し，PCA は共分散行列の固有値分解 $X^\top X = V \Lambda V^\top$ を通じて，分散最大方向への射影を求める \[Ref. 4\]：

```python
sc.tl.pca(adata, svd_solver='arpack')
sc.pl.pca_variance_ratio(adata, log=True)  # elbow plot で PC 数を決定
sc.pl.pca(adata, color=['CST3', 'NKG7', 'PPBP'])
```

通常上位 10〜50 PC を使用する（スクリープロットで決定） \[Ref. 3\]．PCA は線形であるため，scRNA-seq の複雑な非線形構造の可視化には限界がある．

### 近傍グラフの構築

UMAP およびクラスタリングの前に， $k$-最近傍グラフを構築する \[Ref. 1\]：

```python
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
```

### UMAP：位相幾何学的非線形次元削減

UMAP（Uniform Manifold Approximation and Projection; McInnes et al., 2018）は，データが多様体上に分布するという仮定のもと，リーマン幾何学と代数的位相幾何学に基づく理論的枠組みからファジー単体複体（fuzzy simplicial complex）によって高次元の位相的構造をモデル化し，それを低次元で再現する \[Ref. 10\]：

```python
sc.tl.umap(adata)
sc.pl.umap(adata, color=['CST3', 'NKG7', 'PPBP'])
```

アルゴリズムの要点：

1. 高次元空間での $k$-最近傍グラフ構築
2. ファジー重みの計算（指数関数的減衰によるメンバーシップ強度）
3. クロスエントロピーを損失関数として低次元埋め込みを最適化

主要パラメータ \[Ref. 10\]\[Ref. 3\]：

| パラメータ | 効果 | デフォルト |
|---|---|---|
| `n_neighbors` | 局所 vs 大域構造のトレードオフ（小さいと局所的） | 15 |
| `min_dist` | 埋め込み空間での点間最小距離（小さいほど密集） | 0.5 |
| `random_state` | 再現性のための乱数シード | — |

> **数学者への注意**：UMAP の埋め込みは確率的アルゴリズムであり，実行ごとに結果が変わる（`random_state` で固定可能）．また，UMAP は位相的構造を保存するが，ユークリッド的な意味での距離は保存しない．クラスター間の見かけ上の近さや遠さを過度に解釈しないこと \[Ref. 3\]．

---

<a id="module4"></a>
## モジュール 4：クラスタリングとグラフ分割

### Leiden アルゴリズム

scRNA-seq では，まず細胞間の $k$-NN グラフを構築し（モジュール 3 で実施済み），次にグラフ分割によるコミュニティ検出を行う \[Ref. 4\]．**Leiden アルゴリズム**（Traag et al., 2019）は Louvain の改良版で，コミュニティの連結性を保証しつつモジュラリティを最大化する \[Ref. 11\]：

$$Q = \frac{1}{2m} \sum_{ij} \left[ A_{ij} - \frac{k_i k_j}{2m} \right] \delta(c_i, c_j)$$

ここで $A_{ij}$ は隣接行列， $k_i$ は頂点 $i$ の次数， $m$ は総辺数， $\delta(c_i, c_j)$ は $i, j$ が同クラスターのとき 1，そうでないとき 0 となる \[Ref. 11\]：

```python
# igraph 実装（将来の既定）を明示する．n_iterations=2 は Leiden の反復回数，igraph 実装は無向グラフのみ対応
sc.tl.leiden(adata, resolution=0.5, flavor='igraph', n_iterations=2, directed=False)
sc.pl.umap(adata, color=['leiden'], legend_loc='on data')
```

`resolution` パラメータの選択 \[Ref. 3\]\[Ref. 4\]：

- 大きい値 → 細かいクラスター（細かい細胞サブタイプを解像）
- 小さい値 → 粗いクラスター（大きな細胞型ごとの分類）

> **数理的注意**：モジュラリティ最適化は NP 困難問題であり，Leiden アルゴリズムは局所最適解を与える近似アルゴリズムである．`resolution` パラメータの選択に一意の正解はなく，生物学的意味との照合が必要である \[Ref. 11\]．

---

<a id="module5"></a>
## モジュール 5：細胞型アノテーション — マーカー遺伝子と仮説検定

### マーカー遺伝子の探索

各クラスターに特異的に高発現する遺伝子（マーカー遺伝子）を統計的に同定する \[Ref. 4\]\[Ref. 5\]：

```python
sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)
```

Wilcoxon 順位和検定を用いてクラスター内 vs クラスター外の発現量を比較する．これはノンパラメトリック検定であり，正規性の仮定が不要な点が scRNA-seq の歪んだ分布に適している \[Ref. 4\]．

### 既知マーカー遺伝子によるアノテーション

PBMC の代表的なマーカー遺伝子 \[Ref. 4\]\[Ref. 5\]：

| 細胞型 | マーカー遺伝子 | 表面マーカー（参考） |
|---|---|---|
| T 細胞（共通） | CD3D | CD3+ |
| CD4+ T 細胞 | IL7R, CCR7 | CD3+, CD4+ |
| CD8+ T 細胞 | CD8A | CD3+, CD8+ |
| B 細胞 | MS4A1 (CD20), CD79A | CD19+ |
| NK 細胞 | GNLY, NKG7 | CD56+, CD3- |
| 単球（古典的） | CD14, LYZ | CD14+ |
| 単球（非古典的） | FCGR3A, MS4A7 | CD16+ |
| 樹状細胞 | FCER1A, CST3 | — |
| 血小板 | PPBP | — |

```python
marker_genes = ['CD3D', 'IL7R', 'CD8A', 'MS4A1', 'CD79A', 'GNLY', 'NKG7',
                'CD14', 'LYZ', 'FCGR3A', 'MS4A7', 'FCER1A', 'PPBP']
sc.pl.umap(adata, color=marker_genes)
sc.pl.dotplot(adata, marker_genes, groupby='leiden')

# クラスターごとの発現細胞割合（dotplot の点の大きさに相当）
expr = sc.get.obs_df(adata, keys=['leiden'] + marker_genes, use_raw=True)
(expr.set_index('leiden') > 0).groupby(level=0, observed=True).mean().round(2).T
```

PBMC 3k（モジュール 4 の設定）での発現細胞割合の例：

| Leiden | 細胞数 | 主なマーカー（発現細胞割合） | 細胞型 |
|---|---|---|---|
| 0 | 1192 | CD3D 0.88，IL7R 0.66，CD8A 0.12 | CD4 T |
| 1 | 341 | MS4A1 0.86，CD79A 0.94 | B |
| 2 | 641 | CD14 0.54，LYZ 1.00，FCGR3A 0.34，MS4A7 0.40 | Monocyte（CD14+ と FCGR3A+ の混在） |
| 3 | 417 | NKG7 1.00，GNLY 0.57，CD3D 0.59，CD8A 0.35 | CD8 T / NK |
| 4 | 34 | FCER1A 0.82，CST3 1.00 | Dendritic |
| 5 | 13 | PPBP 1.00 | Platelet |

```python
# アノテーションの割り当て（resolution=0.5 で得られた 6 クラスターに対応）
# 根拠は上の dotplot と発現細胞割合の表
cell_type_map = {
    '0': 'CD4 T',       # CD3D, IL7R
    '1': 'B',           # MS4A1, CD79A
    '2': 'Monocyte',    # CD14, LYZ（古典的）と FCGR3A, MS4A7（非古典的）が混在
    '3': 'CD8 T / NK',  # NKG7, GNLY．CD3D・CD8A を発現する細胞（CD8 T）も含む
    '4': 'Dendritic',   # FCER1A, CST3
    '5': 'Platelet',    # PPBP
}
# クラスター構成が変わると対応表は無効になるので確認する
assert set(cell_type_map) == set(adata.obs['leiden'].cat.categories), \
    f"Leiden クラスター {list(adata.obs['leiden'].cat.categories)} が対応表と一致しない"
adata.obs['cell_type'] = adata.obs['leiden'].map(cell_type_map)
sc.pl.umap(adata, color=['cell_type'], legend_loc='on data')
```

> **注意**：モジュール 4 の `resolution=0.5` では Leiden クラスターは 6 個になり，古典的単球（CD14+）と非古典的単球（FCGR3A+），CD8+ T 細胞と NK 細胞がそれぞれ 1 つのクラスターにまとまる．クラスター番号と細胞型の対応は，データ・パラメータ・ソフトウェアの版によって変わる．対応表は既存の例を流用せず，必ずマーカー遺伝子の発現から作り直す．`resolution` を大きくするとこれらのクラスターが分かれるかを確かめるのはよい演習である．

> **自動化ツール**：手動アノテーションが最も難所である．自動化ツールとして ScType（マーカー遺伝子データベース利用）や SingleR（参照データセットとのスピアマン相関）がある \[Ref. 5\]．

---

<a id="module6"></a>
## モジュール 6：差次発現解析（DEG） — 負の二項分布モデル

### なぜ負の二項分布か

scRNA-seq の UMI カウントデータは過分散（分散 > 平均）であり，ポアソン分布ではフィットできない \[Ref. 4\]．遺伝子 $j$ のカウント $X$ は負の二項分布でよく近似される \[Ref. 8\]\[Ref. 9\]：

$$P(X = k) = \binom{k + r - 1}{k} p^k (1-p)^r, \quad k = 0, 1, 2, \ldots$$

平均 $\mu$ と分散 $\sigma^2 = \mu + \mu^2/r$ の関係が二次的であることが，過分散を表現する上で重要である．

### 典型的な比較：細胞型間

```python
# 特定の 2 グループ間の DEG（Wilcoxon 検定）
sc.tl.rank_genes_groups(
    adata, groupby='cell_type',
    groups=['Monocyte'], reference='CD4 T',
    method='wilcoxon'
)
sc.pl.rank_genes_groups(adata, groups=['Monocyte'])
sc.pl.rank_genes_groups_violin(adata, groups='Monocyte', n_genes=8)
```

### 擬似バルク法（Pseudobulk）：集約の実習

複数サンプル（被験者）がある場合，同一サンプル・同一細胞型の細胞のカウントを足し合わせてからサンプル間で比較する（擬似バルク法）．同一被験者由来の細胞は独立ではないため，細胞を独立な観測として検定すると偽陽性が増える（擬似反復問題）\[Ref. 5\]．

**定式化．** 細胞 $i = 1, \dots, n$ のサンプルを $s(i)$，細胞型を $c(i)$ とし，組 $g = (s, c)$ をグループとする．グループ数を $k$ とし，所属行列 $A \in \{0, 1\}^{k \times n}$ を $A_{gi} = 1 \iff (s(i), c(i)) = g$ で定める．生カウント行列 $X \in \mathbb{Z}_{\geq 0}^{n \times m}$ に対し，擬似バルク行列は

$$Y = A X \in \mathbb{Z}_{\geq 0}^{k \times m}, \qquad Y_{gj} = \sum_{i \,:\, (s(i), c(i)) = g} X_{ij}$$

である． $A$ の各列はちょうど 1 つの成分が $1$ なので， $\mathbf{1}_k^\top Y = \mathbf{1}_n^\top X$（遺伝子ごとの総カウントは保存される）．

**注意．**

- 集約するのは **生カウント** である．この時点の `adata.X` は正規化・対数変換・HVG 抽出・スケーリング済みであり，足し合わせても意味がない．そこで生カウントを読み直す．
- PBMC 3k は **単一ドナー** 由来である．ここでは細胞を無作為に $K = 3$ 個の仮サンプルに割り当てるが，これは同一個体の技術的な分割にすぎず，生物学的反復ではない．したがって $Y$ を用いた差次発現検定は行わない（得られる $p$ 値は生物学的推論に使えない）．本来はドナーごとにサンプルを構成する．

```python
import numpy as np
import pandas as pd
import scipy.sparse as sp

# 1. 生カウント行列を読み直し，QC 後の細胞に揃える
#    （この時点の adata.X は正規化・対数変換・HVG 抽出・スケーリング済みで，足し合わせても意味がない）
counts = sc.read_10x_mtx('data/filtered_gene_bc_matrices/hg19/', var_names='gene_symbols')
counts.var_names_make_unique()
counts = counts[adata.obs_names, :].copy()
counts.obs['cell_type'] = adata.obs['cell_type'].values

# 2. 細胞を無作為に K 個の仮サンプルに割り当てる（単一ドナーなので生物学的反復ではない）
K = 3
rng = np.random.default_rng(0)
counts.obs['sample'] = pd.Categorical([f'pseudo{s}' for s in rng.integers(K, size=counts.n_obs)])
print(pd.crosstab(counts.obs['sample'], counts.obs['cell_type']))

# 3. 所属行列 A ∈ {0,1}^{k×n}：A[g, i] = 1 ⇔ 細胞 i がグループ g = (sample, cell_type) に属する
keys = counts.obs['sample'].astype(str) + '|' + counts.obs['cell_type'].astype(str)
group = pd.Categorical(keys)
n = counts.n_obs
k = len(group.categories)
A = sp.csr_matrix((np.ones(n), (group.codes, np.arange(n))), shape=(k, n))

# 4. 擬似バルク行列 Y = A X（k × m）
X = counts.X  # n × m の生カウント（疎行列）
Y = (A @ X).toarray().astype(np.int64)
print(f'X: {X.shape}, A: {A.shape}, Y: {Y.shape}')

pseudobulk = pd.DataFrame(
    Y,
    index=pd.MultiIndex.from_tuples([tuple(g.split('|')) for g in group.categories],
                                    names=['sample', 'cell_type']),
    columns=counts.var_names,
)

# 5. 検算
assert Y.sum() == X.sum()                                   # 総カウントは保存される
n_cells = counts.obs.groupby(['sample', 'cell_type'], observed=True).size()
n_cells.index = n_cells.index.map(lambda t: tuple(map(str, t)))
assert (np.asarray(A.sum(axis=1)).ravel() == n_cells.loc[pseudobulk.index].values).all()  # A の行和 = 各グループの細胞数

agg = sc.get.aggregate(counts, by=['sample', 'cell_type'], func='sum')  # Scanpy 組み込み関数と一致するか
agg_sum = agg.layers['sum']
agg_df = pd.DataFrame(agg_sum.toarray() if sp.issparse(agg_sum) else np.asarray(agg_sum),
                      index=pd.MultiIndex.from_frame(agg.obs[['sample', 'cell_type']].astype(str)),
                      columns=agg.var_names)
assert (agg_df.loc[pseudobulk.index].values == Y).all()
print('検算 OK')

print(f'ゼロ成分の割合: X = {1 - X.nnz / np.prod(X.shape):.3f}, Y = {(Y == 0).mean():.3f}')
pseudobulk[['LYZ', 'CD14', 'MS4A1', 'CD79A', 'IL7R', 'NKG7']]
```

多数のドナーを含む実データで，擬似バルクから負の二項 GLM による検定まで行うのが次の発展版である．

### 発展：複数ドナーでの擬似バルク差次発現解析

ノートブック：[`course/module6-multidonor.ipynb`](../course/module6-multidonor.ipynb)．データは Kang et al. の PBMC（ループス患者 8 名，各ドナーの未処置 `ctrl` と IFN-β 刺激 `stim`）\[Ref. 15\]．配布元の細胞型注釈を用い，CD14+ 単球を解析する．

**モデル．** 擬似バルクサンプル $g$ と遺伝子 $j$ について

$$Y_{gj} \sim \mathrm{NB}(\mu_{gj}, \alpha_j), \qquad \mathrm{Var}(Y_{gj}) = \mu_{gj} + \alpha_j \mu_{gj}^2, \qquad \log \mu_{gj} = \log s_g + \mathbf{x}_g^\top \boldsymbol{\beta}_j$$

とおく（DESeq2 \[Ref. 16\] とその Python 実装 PyDESeq2 \[Ref. 17\] ）．計画式は対応のある設計 `~ replicate + label` とし，条件の係数について $H_0 : \beta_{j,\mathrm{stim}} = 0$ を Wald 検定する．サイズ因子は median-of-ratios

$$s_g = \operatorname*{median}_{j : K_j > 0} \frac{Y_{gj}}{K_j}, \qquad K_j = \Bigl(\prod_{g} Y_{gj}\Bigr)^{1/G}$$

で推定し，分散は遺伝子ごとの推定値をトレンドへ MAP 推定で縮小する．多重検定は Benjamini–Hochberg 法で補正する．ノートブックではトレンドを定数（`fit_type='mean'`）とする．パラメトリックなトレンドの当てはめの収束が遺伝子フィルタの閾値に左右されたためである．

**結果（ $\mathrm{padj} < 0.05$ ， 3709 遺伝子）．**

| 方法 | `stim` 対 `ctrl` で有意な遺伝子数 |
|---|---:|
| 擬似バルク `~ replicate + label` | 2418 |
| 擬似バルク `~ label` | 2352 |
| 細胞単位の Wilcoxon | 1995 |

処置の効果が強いこの比較では，細胞単位の検定の方が有意な遺伝子が多いわけではない．

**擬似反復の実験．** `ctrl` の細胞だけを使い，8 名のドナーを 4 対 4 に分ける $\binom{8}{4} / 2 = 35$ 通りすべてで，群 A，B を比較した．群に効果はないので，有意な遺伝子は偽陽性である．

| 方法 | 偽陽性 0 個の分け方 | 偽陽性数の中央値 | 偽陽性数の最大値 |
|---|---:|---:|---:|
| 細胞単位の Wilcoxon | 0 / 35 | 139 | 271 |
| 擬似バルク DESeq2 | 27 / 35 | 0 | 24 |

細胞単位の検定はドナー間の差を群の差と取り違える \[Ref. 18\]\[Ref. 19\]．擬似バルクは偽陽性を大幅に減らすが，なくすわけではない．全帰無のもとで BH 法の FDR は「1 個以上棄却する確率」に等しく，独立または PRDS の仮定のもとで $0.05$ 以下のはずだが，観測は $8 / 35 \approx 0.23$ だった．分け方どうしは独立でないので正確な率とは言えず，原因（小標本での Wald 検定の正規近似など）は確かめていない．

---

<a id="module7"></a>
## モジュール 7：発展 — 軌跡推定と擬似時間

本モジュールは PBMC データよりも，骨髄造血や神経分化のような連続的変化を示す発達系データでより意義があるが，概念学習として PBMC でも実施可能である．

### PAGA（Partition-based Graph Abstraction）

細胞状態の連続的変化（分化過程など）を推定する手法（Wolf et al., 2019）\[Ref. 12\]．クラスター間のコネクティビティをグラフの辺の重みとして定量化し，細胞状態空間の位相的構造を可視化する：

```python
sc.tl.paga(adata, groups='leiden')
sc.pl.paga(adata, plot=True, threshold=0.03)
```

### 拡散擬似時間（Diffusion Pseudotime, DPT）

拡散写像（Diffusion Maps; Coifman & Lafon, 2006）\[Ref. 13\] に基づき，参照細胞（根細胞）からの擬似的な「発生時間」を計算する \[Ref. 4\]：

```python
import numpy as np

# 根細胞（幹細胞に近い細胞）を指定
adata.uns['iroot'] = np.flatnonzero(adata.obs['leiden'] == '0')[0]
sc.tl.diffmap(adata)  # 拡散写像（DPT の前提）を先に計算する
sc.tl.dpt(adata)
sc.pl.umap(adata, color=['dpt_pseudotime'], color_map='viridis')
```

### RNA Velocity（参考）

スプライシング前後の mRNA 比率から細胞の発現状態の「速度ベクトル」を推定する手法が **scVelo**（Bergen et al., 2020）である \[Ref. 14\]．これは常微分方程式系のパラメータ推定問題として定式化できる：

$$\frac{du}{dt} = \alpha(t) - \beta u, \qquad \frac{ds}{dt} = \beta u - \gamma s$$

ここで $u$ は未スプライシング mRNA， $s$ はスプライシング済み mRNA， $\alpha$ は転写率， $\beta$ はスプライシング率， $\gamma$ は分解率である \[Ref. 14\]：

> **⚠ データの準備**：scVelo は unspliced（スプライシング前）と spliced（スプライシング後）の両カウントを必要とする．10x Genomics の `filtered_gene_bc_matrices` にはこの情報が含まれないため，**velocyto** または **STARsolo** を用いて別途 loom / h5ad ファイルを生成する必要がある．PBMC 3k でそのまま実行することはできない．

本コースでは，scVelo が提供する組み込みサンプルデータ（膵臓β細胞分化）を使用して動作を確認する：

```python
import scvelo as scv

# scVelo サンプルデータ（膵臓β細胞分化）を使用
# unspliced / spliced カウントが含まれた AnnData が返される
adata_velo = scv.datasets.pancreas()

scv.pp.filter_and_normalize(adata_velo, min_shared_counts=20)
# 近傍グラフは Scanpy で先に作る（scVelo による自動計算は非推奨）
sc.pp.pca(adata_velo, n_comps=30)
sc.pp.neighbors(adata_velo, n_neighbors=30, n_pcs=30)
scv.pp.moments(adata_velo, n_pcs=None, n_neighbors=None)

# dynamical モデルの前に splicing dynamics を推定する
scv.tl.recover_dynamics(adata_velo)
scv.tl.velocity(adata_velo, mode='dynamical')
scv.tl.velocity_graph(adata_velo)
scv.pl.velocity_embedding_stream(adata_velo, basis='umap', color='clusters')
```

> **注意（scvelo 0.3.x 以降）**：旧来の `scv.read()` は廃止されており，`scv.load()` を使用すること．また `filter_and_normalize()` の `n_top_genes` 引数も削除されているため使用不可．`mode='dynamical'` を指定する場合は事前に `scv.tl.recover_dynamics()` の実行が必要である．

自前のデータ（velocyto / STARsolo 出力の loom ファイル）を使う場合：

```python
import scvelo as scv

# loom ファイルの読み込み（scv.read は廃止，scv.load を使用）
adata_velo = scv.load('path/to/sample.loom', cache=True)
scv.utils.clean_obs_names(adata_velo)

# モジュール 1〜5 で作成した adata とバーコードを照合してマージ
adata_velo = scv.utils.merge(adata_velo, adata)

scv.pp.filter_and_normalize(adata_velo, min_shared_counts=20)
# 近傍グラフは Scanpy で先に作る（scVelo による自動計算は非推奨）
sc.pp.pca(adata_velo, n_comps=30)
sc.pp.neighbors(adata_velo, n_neighbors=30, n_pcs=30)
scv.pp.moments(adata_velo, n_pcs=None, n_neighbors=None)
scv.tl.recover_dynamics(adata_velo)
scv.tl.velocity(adata_velo, mode='dynamical')
scv.tl.velocity_graph(adata_velo)
scv.pl.velocity_embedding_stream(adata_velo, basis='umap')
```

---

## 数学者向け概念マッピング

scRNA-seq 解析の各ステップと数理的概念の対応を整理する：

| 解析ステップ | 対応する数理概念 | 参照 |
|---|---|---|
| カウント行列 | 非負整数値疎行列 $\mathbb{Z}_{\geq 0}^{m \times n}$ | \[Ref. 1\] |
| QC フィルタリング | 外れ値検出（ロバスト統計） | \[Ref. 4\]\[Ref. 5\] |
| 正規化 | スケーリング変換，対数変換 | \[Ref. 8\]\[Ref. 9\] |
| HVG 選択 | 分散の正則化推定，特徴選択 | \[Ref. 4\] |
| PCA | 特異値分解（SVD），主成分分析 | \[Ref. 4\] |
| $k$-NN グラフ | 計量空間上のグラフ構築 | \[Ref. 10\] |
| UMAP | 位相的データ解析，多様体学習 | \[Ref. 10\] |
| Leiden クラスタリング | グラフのモジュラリティ最適化 | \[Ref. 11\] |
| Wilcoxon 検定 | ノンパラメトリック仮説検定 | \[Ref. 4\] |
| 負の二項分布モデル | 過分散カウントデータの確率モデル | \[Ref. 8\] |
| 拡散擬似時間 | マルコフ連鎖，拡散過程 | \[Ref. 13\] |
| RNA velocity | 常微分方程式系，パラメータ推定 | \[Ref. 14\] |

---

## 推奨学習リソースまとめ

| リソース | 種別 | URL |
|---|---|---|
| Scanpy PBMC 3k 公式チュートリアル（legacy workflow，Scanpy 1.12.x） | Jupyter Notebook | https://scanpy.readthedocs.io/en/1.12.x/tutorials/basics/clustering-2017.html |
| Single-cell best practices（オンライン教科書） | Web/Jupyter | https://www.sc-best-practices.org |
| theislab/single-cell-tutorial（GitHub） | コード | https://github.com/theislab/single-cell-tutorial |
| Galaxy Training: Clustering 3K PBMCs with Scanpy | ブラウザ実行可 | https://training.galaxyproject.org/training-material/topics/single-cell/tutorials/scrna-scanpy-pbmc3k/tutorial.html |
| Luecken & Theis (2019) 論文 | 査読済み論文 | https://doi.org/10.15252/msb.20188746 |
| Heumos et al. (2023) 論文 | 査読済み論文 | https://www.nature.com/articles/s41576-023-00586-w |

---

<a id="references"></a>
## 参考文献

**\[Ref. 1\]** Wolf, F.A., Angerer, P., Theis, F.J. (2018). SCANPY: large-scale single-cell gene expression data analysis. *Genome Biology*, 19, 15. https://doi.org/10.1186/s13059-017-1382-0

**\[Ref. 2\]** Preprocessing and clustering 3k PBMCs (legacy workflow). Scanpy Documentation (version 1.12.x). https://scanpy.readthedocs.io/en/1.12.x/tutorials/basics/clustering-2017.html

**\[Ref. 3\]** Single-cell best practices. Theis Lab, Helmholtz Munich. Online book (continuously updated). https://www.sc-best-practices.org

**\[Ref. 4\]** Luecken, M.D., Theis, F.J. (2019). Current best practices in single-cell RNA-seq analysis: a tutorial. *Molecular Systems Biology*, 15(6), e8746. https://doi.org/10.15252/msb.20188746

**\[Ref. 5\]** Heumos, L., Schaar, A.C., Lance, C., Litinetskaya, A., Drost, F., Zappia, L., ... Theis, F.J. (2023). Best practices for single-cell analysis across modalities. *Nature Reviews Genetics*, 24, 550–572. https://doi.org/10.1038/s41576-023-00586-w

**\[Ref. 6\]** Young, M.D., Behjati, S. (2020). SoupX removes ambient RNA contamination from droplet-based single-cell RNA sequencing data. *GigaScience*, 9(12), giaa151. https://doi.org/10.1093/gigascience/giaa151

**\[Ref. 7\]** Germain, P.-L., Lun, A., Garcia Meixide, C., Macnair, W., Robinson, M.D. (2022). Doublet identification in single-cell sequencing data using scDblFinder. *F1000Research*, 10, 979. https://doi.org/10.12688/f1000research.73600.2

**\[Ref. 8\]** Hafemeister, C., Satija, R. (2019). Normalization and variance stabilization of single-cell RNA-seq data using regularized negative binomial regression. *Genome Biology*, 20, 296. https://doi.org/10.1186/s13059-019-1874-1

**\[Ref. 9\]** Lun, A.T.L., Bach, K., Marioni, J.C. (2016). Pooling across cells to normalize single-cell RNA sequencing data with many zero counts. *Genome Biology*, 17, 75. https://doi.org/10.1186/s13059-016-0947-7

**\[Ref. 10\]** McInnes, L., Healy, J., Melville, J. (2018). UMAP: Uniform Manifold Approximation and Projection for Dimension Reduction. *arXiv*:1802.03426. https://arxiv.org/abs/1802.03426

**\[Ref. 11\]** Traag, V.A., Waltman, L., van Eck, N.J. (2019). From Louvain to Leiden: guaranteeing well-connected communities. *Scientific Reports*, 9, 5233. https://doi.org/10.1038/s41598-019-41695-z

**\[Ref. 12\]** Wolf, F.A., Hamey, F.K., Plass, M., Solana, J., Dahlin, J.S., Göttgens, B., ... Theis, F.J. (2019). PAGA: graph abstraction reconciles clustering with trajectory inference through a topology preserving map of single cells. *Genome Biology*, 20, 59. https://doi.org/10.1186/s13059-019-1663-x

**\[Ref. 13\]** Coifman, R.R., Lafon, S. (2006). Diffusion maps. *Applied and Computational Harmonic Analysis*, 21(1), 5–30. https://doi.org/10.1016/j.acha.2006.04.006

**\[Ref. 14\]** Bergen, V., Lange, M., Peidli, S., Wolf, F.A., Theis, F.J. (2020). Generalizing RNA velocity to transient cell states through dynamical modeling. *Nature Biotechnology*, 38, 1408–1414. https://doi.org/10.1038/s41587-020-0591-3

**\[Ref. 15\]** Kang, H.M., Subramaniam, M., Targ, S., Nguyen, M., Maliskova, L., ... Ye, C.J. (2018). Multiplexed droplet single-cell RNA-sequencing using natural genetic variation. *Nature Biotechnology*, 36(1), 89–94. https://doi.org/10.1038/nbt.4042

**\[Ref. 16\]** Love, M.I., Huber, W., Anders, S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology*, 15, 550. https://doi.org/10.1186/s13059-014-0550-8

**\[Ref. 17\]** Muzellec, B., Teleńczuk, M., Cabeli, V., Andreux, M. (2023). PyDESeq2: a python package for bulk RNA-seq differential expression analysis. *Bioinformatics*, 39(9), btad547. https://doi.org/10.1093/bioinformatics/btad547

**\[Ref. 18\]** Squair, J.W., Gautier, M., Kathe, C., Anderson, M.A., ... Courtine, G. (2021). Confronting false discoveries in single-cell differential expression. *Nature Communications*, 12, 5692. https://doi.org/10.1038/s41467-021-25960-2

**\[Ref. 19\]** Zimmerman, K.D., Espeland, M.A., Langefeld, C.D. (2021). A practical solution to pseudoreplication bias in single-cell studies. *Nature Communications*, 12, 738. https://doi.org/10.1038/s41467-021-21038-1

