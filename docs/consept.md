## 最終目標
- カセットエキソンを抽出し、target AidなどのBEでエキソンスキップが可能かどうか出力する
- 遺伝子名を打つとトランスクリプトが出力される
- そのトランスクリプトを消せるsgRNAを出力する

## 考察
スプライスバリアントがある遺伝子はそもそも何パーセント？
バリアント数ごとに、何個壊せるかを出力する。
5'の方は転写開始点なので、プロモーターエンハンサーの位置が遠くなったりすると考えると、UTRだけが違っても意味はありそう。  
A3SSとA5SSが起こる場合は、コーディングシーケンス中にスプライスアクセプター・ドナーが存在することになるので、そこを編集すると残したいバリアントも同時に消えてしまうことになる  


## 小タスクに分解する
- [ ] refflatのデータ構造を把握する  
geneName name(=transcript id) chrom strand  txStart txEnd  cdsStart cdsEnd exonCount exonStarts exonEnds  

- [ ] 重複列がないことを確認する  

- [ ] エキソンのスタートエンドを　[[start1,end1],[start2,end2],....]　構造に変換する    
  
  modify_refflat()  

- [ ] エキソンの種類を判定する関数の作成    
  入力: 調べたいexonのstart,endペア,ある遺伝子のすべてのトランスクリプトが持つエキソンのstart end のペア(target_exon, all_transcripts, fuzzy_tolerance=5)    

　判定の基準: ある(start:end) ペアに対して、  

  - startはすべてのtranscriptに存在し、endもすべてのtranscriptに存在する → "constitutive"  
  - 他のtranscriptの（start: end）の組に対して、startは一致しているが、endだけ異なる場合がある  → "a5ss" (次のエキソンから見た時の5'側のスプライスサイト変化)  
  - 他のtranscriptの（start: end）の組に対して、endは一致しているが、startは一致しない → "a3ss" (前のエキソンから見た時の3'側のスプライスサイト変化)
  - startはすべてのtranscriptに存在せず、endもすべてのtransctiptに存在しない。しかし、start:endの組が2つ以上のtranscriptに存在する→ "cassette"  
  - startはすべてのtranscriptに存在せず、endもすべてのtransctiptに存在しない。そして、start:endの組が1つのtranscriptにしか存在しない→ "unique"  
  - start:endが他のtranscriptと一致しないが、他のエキソンと1塩基以上のオーバーラップが生じているもの→ "overlap" 

  出力: [constitutive, cassette, unique....]
  
  classify_exon_type()

- [ ] 遺伝子ごとにrefflatをgroup化して、それぞれのエキソンに対してclassify_exon_type()を実行し、"exon_type"列に結果を格納する

--- データ解析 ---  

- [ ] 各エキソンカテゴリごとの分布などを解析する
- [ ] それらの抽出したエキソン長は3の倍数であることを確認する（3の倍数でないとエキソンスキップによってフレームシフトが起こり得るため、カセットエキソンになり得ない？）  
- [ ] uniqueのうち、エキソン長が3の倍数のエキソンのみをカセットエキソンとして扱う

- [ ] 取得したカセットエキソンの位置情報をマウスゲノムFASTAにマッピングする  
- [ ] エキソン開始点・終着点の+-25bp程度を取得する  

--- 未定 ----  
- [ ] TargetAID (SpCas9, NGG)の編集ウィンドウに入るsgRNAをデザインする