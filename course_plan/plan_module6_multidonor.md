# 計画：モジュール 6 発展版 — 複数ドナーデータによる擬似バルク差次発現解析

> **状態**：計画（未着手）．記録日 2026-11-08．
> **位置づけ**：モジュール 6 の擬似バルク節（PBMC 3k の無作為分割による「集約の実習」）の続き．

## 背景と目的

PBMC 3k は単一ドナー由来のため，擬似バルク法の本来の目的（被験者間のばらつきを考慮した条件間比較）を体験できない．現行のモジュール 6 では，無作為に作った仮サンプルで集約 $Y = AX$ を計算するところまでにとどめ，検定は行っていない．

発展版では，**生物学的反復（複数ドナー）と条件（処置／対照）を持つ公開データ** を用いて，

1. ドナー × 細胞型 × 条件 ごとの擬似バルク行列を作り，
2. 負の二項一般化線形モデル（DESeq2 の Python 実装 PyDESeq2）で条件間の差次発現を検定し，
3. 細胞単位の Wilcoxon 検定（擬似反復）との結果の違いを比較する

ことで，モジュール 6 の主題である **負の二項分布モデル** と **擬似反復問題** を実データで扱う．

## データ候補

**Kang et al. (2018) の IFN-β 刺激 PBMC データ**

- 文献：Kang, H.M., Subramaniam, M., Targ, S., et al. (2018). Multiplexed droplet single-cell RNA-sequencing using natural genetic variation. *Nature Biotechnology*, 36, 89–94. https://doi.org/10.1038/nbt.4042 （Crossref で書誌確認済み．オンライン公開は 2017-12-11．Author Correction: https://doi.org/10.1038/s41587-020-0715-9）
- 内容：ループス患者 8 名の PBMC．各ドナーについて未処置（`ctrl`）と IFN-β 刺激（`stim`）の 2 条件（sc-best-practices の差次発現の章で確認済み）．
- `obs` の主な列：`label`（条件），`replicate`（患者 ID），`cell_type`（sc-best-practices の章の出力で確認済み）．
- 対応のある設計：同一ドナーの `ctrl` / `stim` が対になる．

## 参考にする一次資料

- Single-cell best practices，「Differential gene expression analysis」の章：https://www.sc-best-practices.org/conditions/differential-gene-expression/ （Kang データ，擬似バルク，PyDESeq2 を使用．pertpy チュートリアルを使えば edgeR 版にも書き換えられると記載）
- PyDESeq2：https://pypi.org/project/pydeseq2/ （2026-11-08 時点で 0.5.4，Python >= 3.11）
- pertpy の `pertpy.data.kang_2018`：https://pertpy.readthedocs.io/en/latest/api/data/pertpy.data.kang_2018.html

## 手順案（ノートブックの新しい節）

1. **データ取得**：Kang データを AnnData として読み込み，生カウントであることを確認する．
2. **擬似バルク化**：`(replicate, label, cell_type)` ごとに和をとる．既存節の所属行列 $A$ による定式化，または `sc.get.aggregate` を用いる．細胞数の少ないグループは除外する（閾値は資料に合わせて決める）．
3. **モデル**：1 つの細胞型（例：CD14 単球）に絞り，遺伝子 $j$・擬似バルクサンプル $g$ について
   $$Y_{gj} \sim \mathrm{NB}(\mu_{gj}, \alpha_j), \qquad \log \mu_{gj} = \log s_g + \mathbf{x}_g^\top \boldsymbol{\beta}_j$$
   とおく．計画行列 $\mathbf{x}_g$ にはドナー効果と条件効果を入れる（例：`~ replicate + label`）．$s_g$ はサイズ因子，$\alpha_j$ は分散パラメータ．DESeq2 のサイズ因子推定，分散の縮小推定，Wald 検定を数式で解説する．
4. **検定**：PyDESeq2 で `stim` 対 `ctrl` を検定し，調整済み $p$ 値と log fold change を得る．
5. **比較**：同じ細胞型について細胞単位の Wilcoxon 検定（`sc.tl.rank_genes_groups`）を行い，有意遺伝子数と $p$ 値の分布を比べる．擬似反復で偽陽性が増える様子を確認する．
6. **（任意）**：decoupler によるパスウェイ活性解析（sc-best-practices の続きの章）．

## 決めておくこと・未確認事項

- [ ] **データの取得方法**：pertpy 1.5.0（2026-11-08 時点の最新）は Python >= 3.13 を要求する．現在の `course/.venv` は Python 3.12 なので，次のどれかを選ぶ．
  - (a) 環境を Python 3.13 に上げて pertpy を入れる．
  - (b) pertpy の古い版で Python 3.12 に対応するものを使う（対応版は未確認）．
  - (c) pertpy を使わず，データファイル（h5ad）を直接ダウンロードする（配布元 URL は未確認）．
- [ ] PyDESeq2 0.5.x の API（`DeseqDataSet` / `DeseqStats` の引数名，計画式の書き方）を公式ドキュメントで確認する．
- [ ] 擬似バルクにするときの最小細胞数の閾値を，sc-best-practices の章に合わせて決める．
- [ ] 追加するパッケージを `requirements.txt` に反映する（`uv pip freeze > requirements.txt`）．
- [ ] 参考文献（Kang et al. 2018，PyDESeq2，Squair et al. 2021 など擬似反復の文献）を `scrna_hands_on_course.md` の参考文献リストに追加する．書誌情報は Crossref 等で確認する．
- [ ] データのダウンロードサイズと実行時間を見積もる．
