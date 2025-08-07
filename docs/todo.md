## 最終目標
- Skipped exon (SE) を抽出し、target AidなどのBEでエキソンスキップが可能かどうか出力する
- 遺伝子名を入力するとトランスクリプトが出力される
- 各skippedまたはuniqueエキソンのSA・SDを編集すsgRNAを設計し、表示する

## 考察
スプライスバリアントがある遺伝子はそもそも何パーセント？
バリアント数ごとに、何個壊せるかを出力する。
5'の方は転写開始点なので、プロモーターエンハンサーの位置が遠くなったりすると考えると、UTRだけが違っても意味はありそう。  
A3SSとA5SSが起こる場合は、コーディングシーケンス中にスプライスアクセプター・ドナーが存在することになるので、そこを編集すると残したいバリアントも同時に消えてしまうことになる  


## データ構造の把握、前処理
- [x] refflatのデータ構造を把握する  
geneName name(=transcript id) chrom strand  txStart txEnd  cdsStart cdsEnd exonCount exonStarts exonEnds  

- [x] 重複列がないことを確認する。トランスクリプト名に重複が存在するため、重複列を完全に削除する  

- [x] エキソンのスタートエンドを　[[start1,end1],[start2,end2],....]　構造に変換し、エキソン長を計算する   
   

## スプライシングイベントに応じたエキソンへのアノテーション付加
- [x] エキソンの種類を判定する関数の作成    
  入力: 調べたいexonのstart,endペア,ある遺伝子のすべてのトランスクリプトが持つエキソンのstart end のペア(target_exon, all_transcripts）

　判定の基準: ある(start:end) ペアに対して、  

  - "constitutive" : startはすべてのtranscriptに存在し、endもすべてのtranscriptに存在する  
  - "a5ss": 他のtranscriptの（start: end）の組に対して、startは一致しているが、endだけ異なる場合がある
    - "a5ss-short": startを共有するエキソンの集団で、自分よりendが大きいエキソンが存在する
    - "a3ss-long": startを共有するエキソンの集団で、自分が最も大きいendを持つ
  - "a3ss": 他のtranscriptの（start: end）の組に対して、endは一致しているが、startは一致しない 
    - "a5ss-short": endを共有するエキソンの集団で、自分よりstartが小さいエキソンが存在する
    - "a5ss-long": endを共有するエキソンの集団で、自分が最も小さいstartを持つ
  - "skipped exon": startはすべてのtranscriptに存在せず、endもすべてのtransctiptに存在しない。しかし、start:endの組が2つ以上のtranscriptに存在する  
  - "unique": startはすべてのtranscriptに存在せず、endもすべてのtransctiptに存在しない。そして、start:endの組が1つのtranscriptにしか存在しない
  - "overlap":start: endが他のtranscriptと一致しないが、他のエキソンと1塩基以上のオーバーラップが生じているもの 
  - "intron_retention": あるエキソンに対してa3ssとa5ssに相当するエキソンが両方存在するもの

  出力: [constitutive, skipped_exon, unique....]

- [x] 遺伝子ごとにrefflatをgroup化して、それぞれのエキソンに対して```classify_exon_type()```を実行し、```"exon_type"```列に結果を格納する  

## フレームやコーディングなどのアノテーションを付加する
- [x] 常染色体または性染色体だけにマッピングされている遺伝子を残し、そうでないものは削除する
- [x] エキソン長が3の倍数かどうかのアノテーションを"flame" 列に追加する（out-flameまたはin-flame)
- [x] その遺伝子がタンパク質をコードするかの情報を"cording"列に追加する（cordingまたはnoncording）
- [x] エキソンの位置を"exon_position"列に追加する(first, internal, last) →　ある時はfirst, ある時はinternal のような場合は存在するのか？

## データ解析
- [x] 全遺伝子に対してエキソンの各アノテーションがどれくらい存在するかを調べる
- [x] Skipped exonが存在し得ないトランスクリプト（スプライシングバリアントが存在しない、またはエキソン数が1の遺伝子）を除いて各アノテーションを持つ遺伝子数を調べる
- [x] Skipped or Uniqueが中間エキソンまたは端エキソンのどちらに存在するか調べてみる
- [x] 取得したSA/SD周辺配列がどのくらいcanonical splice siteかを調べる

## 塩基配列取得のための前処理
- [x] Skipped or Unique or A5ss-long or A3ss-long にアノテーションされたエキソンを少なくとも一つ持つトランスクリプトだけをフィルターする
- [x] 1エキソンー1行になるようにrefflatを展開する
- [x] ターゲットのエキソンだけを抽出する（各エキソンに一意であるUUIDをscore列に追加する）（重複を削除する）
- [x] それらのExonStart, ExonEndから+-25bp(トータルで50bp)の範囲を抜き出し、保存する
- [x] dfをBED形式に変換し、保存する

## ゲノム配列の取得
- [x] 取得したカセットエキソンの位置情報をマウスゲノムにマッピングし、塩基配列を取得する(pybedtoolを利用)
- [x] 塩基配列を付加したDFを作成する
- [x] 任意のPAMと編集ウィンドウを持つCBEに対してsgRNAをデザインする
- [ ] 任意のPAMと編集ウィンドウを持つABEに対してsgRNAをデザインする
- [ ] よく使われているBE + target-AIDの条件で関数を実行して結果をdfに付加する

## デザインしたデータの可視化
- [ ] -view モードで遺伝子名をクエリに、検索できるようにする
- [ ] IGVやUCSC genome browser custom track に投げられるフォーマットに変更する 
