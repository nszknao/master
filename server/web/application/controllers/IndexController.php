<?php
require_once ("controllers/Abstract.php");

/**
 * IndexController
 * @author iLoiLohas
 */
class
	IndexController
extends
	ControllerAbstract
{

	public function init()
	{
		$this->_loginit(get_class($this));
	}
	public function preDispatch() {
		/* 各コントローラの共通前処理 */
		parent::_preDispatch();
	}
	/**
	 * ログイン処理
	 * route --> /
	 */
	public function indexAction() {
		$this->_log->debug(__CLASS__ . ":" . __FUNCTION__ . " called:(" . __LINE__ . ")");
	}
	/**
	 * 計算開始
	 * route --> /run
	 */
	public function runAction() {
		$this->_log->debug(__CLASS__ . ":" . __FUNCTION__ . " called:(" . __LINE__ . ")");
		$params	= $this->getAllParams();
		if (count($params) == 0) return;

		switch ($params['simulation']) {
			case 'both':
				$output	= shell_exec(APPLICATION_PATH.'/../../sh/run.sh');
				break;
			case 'simulation':
				$output	= exec(APPLICATION_PATH.'/../../sh/sRun.sh');
				if ($output == NULL) {
					$this->_log->debug("エラーが発生しました。");
					break;
				}
				$this->_log->debug("sRun.shを実行しました。");
				$this->_log->debug(print_r($output, true));
				break;
			case 'analysis':
				$output	= shell_exec(APPLICATION_PATH.'/../../sh/aRun.sh');
				break;
		}
		$this->redirect("/");
		return;
	}
}