<?php
/**
 * PHPの環境設定
 * @author iLoiLohas
 */
set_include_path(
	implode(PATH_SEPARATOR
	,	array(
			APPLICATION_PATH
		,	APPLICATION_PATH . '/../library/log4php'
		,	APPLICATION_PATH . '/../library/Smarty'
		,	APPLICATION_PATH . '/../library/vendor'
		,	get_include_path(),
		)
	)
);