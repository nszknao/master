{include file='common/default_css.tpl' assign=default_css}
{include file='common/default_js.tpl' assign=default_js}

<!DOCTYPE HTML>
<html lang="ja">
<head>
	<meta charset="UTF-8">
	<title></title>
{$default_css}
</head>
<body>
	<div class="container">
		<div class="page-header">
			<h3>計算を実行する</h3>
		</div>
		<form action="/run" method="post" class="form-horizontal">
			<div class="form-group">
				<label class="radio-inline col-sm-3 col-sm-offset-1"><input type="radio" name="simulation">解析＋シミュレーション</label>
				<label class="radio-inline col-sm-3 col-sm-offset-1"><input type="radio" name="simulation">シミュレーション</label>
				<label class="radio-inline col-sm-3 col-sm-offset-1"><input type="radio" name="simulation">解析</label>
			</div>
			<div class="form-group">
				<div class="col-sm-offset-9">
					<button type="submit" class="btn btn-default">実行</button>
				</div>
			</div>
		</form>
	</div>
{$default_js}
</body>
</html>