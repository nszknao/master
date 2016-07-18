<?php
class Common
{
	/**
	 * マスターDBを取得する．
	 * @return Zend_Db_Adapter_Mysqli
	 */
	static public function getMaster() {
		$db	= Zend_Registry::get('db1');
		$db->getConnection();
		$db->setFetchMode(Zend_DB::FETCH_ASSOC);
		return $db;
	}
	/**
	 * スレイブDBを取得する．
	 * @return Zend_Db
	 */
	static public function getSlave() {
		$db	= Zend_Registry::get('db2');
		$db->getConnection();
		$db->setFetchMode(Zend_DB::FETCH_ASSOC);
		return $db;
	}
	/**
	 * ワークセッション名を取得する。
	 */
	public static function  getSessionName() {
		$space	= "login_".@$_SERVER['HTTP_HOST'];
		return $space;
	}
}