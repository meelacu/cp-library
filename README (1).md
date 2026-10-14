# cp-library

競技プログラミング用の C++ ライブラリ。

## 生成AIの利用について / Notice on generative AI

このリポジトリには、生成AIの支援を受けて作成したコードと個人で作成したコードが含まれます。
AtCoder の開催中のコンテスト（ABC / ARC / AGC）でこのリポジトリのコードを使う場合は、
[AtCoder生成AI対策ルール](https://info.atcoder.jp/entry/llm-rules-ja) に従い、使用したファイルの URL を提出コードに記載してください。

This repository contains code written with the help of generative AI.
If you use it during an AtCoder contest (ABC / ARC / AGC), include the URL of the file you used in your submission,
as required by [AtCoder's rules against generative AI](https://info.atcoder.jp/entry/llm-rules-en).

## 構成

```
lib/       ライブラリ本体（1ファイル = 1機能）
snip.py    lib/ から VSCode のスニペットを生成するスクリプト
```

## 使い方

`lib/` のファイルの中身をそのままコピーして使えます。

- 先頭の `// @prefix` などの行はスニペット生成用の情報で、コードとしては不要です。
- 本体中の `$0` や `${1:...}` は VSCode スニペットのプレースホルダです。

VSCode のスニペットとして使う場合は、Python 3.8 以上で次を実行すると、`lib/` から生成したスニペットがユーザースニペットのフォルダに書き出されます。

```
python snip.py build
```

各スニペットの先頭には、そのファイルの GitHub 上の URL がコメントとして入ります（fork した場合は fork 先の URL）。
