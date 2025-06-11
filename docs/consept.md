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
  
  ```modify_refflat()```  

## スプライシングイベントに応じたエキソンへのアノテーション付加
- [x] エキソンの種類を判定する関数の作成    
  入力: 調べたいexonのstart,endペア,ある遺伝子のすべてのトランスクリプトが持つエキソンのstart end のペア(target_exon, all_transcripts）

　判定の基準: ある(start:end) ペアに対して、  

  - startはすべてのtranscriptに存在し、endもすべてのtranscriptに存在する → "constitutive"  
  - 他のtranscriptの（start: end）の組に対して、startは一致しているが、endだけ異なる場合がある  → "a5ss" (次のエキソンから見た時の5'側のスプライスサイト変化)  
  - 他のtranscriptの（start: end）の組に対して、endは一致しているが、startは一致しない → "a3ss" (前のエキソンから見た時の3'側のスプライスサイト変化)
  - startはすべてのtranscriptに存在せず、endもすべてのtransctiptに存在しない。しかし、start:endの組が2つ以上のtranscriptに存在する→ "skipped exon"  
  - startはすべてのtranscriptに存在せず、endもすべてのtransctiptに存在しない。そして、start:endの組が1つのtranscriptにしか存在しない→ "unique"  
  - start:endが他のtranscriptと一致しないが、他のエキソンと1塩基以上のオーバーラップが生じているもの→ "overlap"
  - あるエキソンに対してa3ssとa5ssに相当するエキソンが両方存在するもの→ "split"

  出力: [constitutive, skipped_exon, unique....]
  
  ```classify_exon_type()```

- [x] 遺伝子ごとにrefflatをgroup化して、それぞれのエキソンに対してclassify_exon_type()を実行し、"exon_type"列に結果を格納する  

## フレームやコーディングなどのアノテーションを付加する
- [x] 常染色体または性染色体だけにマッピングされている遺伝子を残し、そうでないものは削除する
- [ ] エキソン長が3の倍数かどうかのアノテーションを"flame" 列に追加する（out-flameまたはin-flame)
- [ ] その遺伝子がタンパク質をコードするかの情報を"cording"列に追加する（cordingまたはnoncording）

## データ解析
- [ ] 全遺伝子に対してエキソンの各アノテーションがどれくらい存在するかを調べる
- [ ] Skipped exonが存在し得ないトランスクリプト（スプライシングバリアントが存在しない、またはエキソン数が1の遺伝子）を除いて各アノテーションを持つ遺伝子数を調べる
- [ ] skipped or uniqueが中間エキソンまたは端エキソンのどちらに存在するか調べて

## ゲノム配列の取得
- [ ] 取得したカセットエキソンの位置情報をマウスゲノムFASTAにマッピングする  
- [ ] エキソン開始点・終着点の+-25bp程度を取得する  


- [ ] TargetAID (SpCas9, NGG)の編集ウィンドウに入るsgRNAをデザインする
