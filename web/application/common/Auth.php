<?php
require_once 'common/Common.php';

class
	Auth
{
	/**
	 * ログインのチェックを行う．
	 */
	public static function loginCheck() {
		$sess	= self::_getAuth();
		if (count($sess) == 0) {
			return false;
		}
		return $sess;
	}
	/**
	 * Auth情報をセッションから取得する．
	 */
	protected function _getAuth() {
		Zend_Session::start();
		$sess	= new Zend_Session_Namespace(Common::getSessionName());
		if (!isset($sess->login)) {
			$sess->login	= array();
		}
		return $sess->login;
	}
	/**
	 * Auth情報をセッションへ設定する．
	 */
	public static function setAuth($user) {
		$sess	= new Zend_Session_Namespace(Common::getSessionName());
		$sess->login	= $user;
		$sess->setExpirationSeconds(600);
	}
	/**
	 * Auth情報をセッションから削除する．
	 */
	public static function deleteAuth() {
		$sess	= new Zend_Session_Namespace(Common::getSessionName());
		$sess->login	= array();
	}
}