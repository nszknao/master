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
		$this->redirect("/");
	}
}