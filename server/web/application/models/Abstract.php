<?php
require_once 'common/Common.php';
require_once 'common/Auth.php';

class ModelAbstract
{
	protected $_log	= NULL;
	protected $_dbs	= array();
	
	/**
	 *
	 * ログ初期化.
	 * @return Logger
	 */
	protected function _loginit($classname) {
		$log		= Zend_Registry::get('logger');
		$this->_log	= $log::getLogger($classname);
		return	$this->_log;
	}
	/**
	 * データベーストランザクションを開始する．
	 * @param Zend_Db_Adapter_Mysqli $db
	 */
	protected function _begin($db) {
		$this->_log->debug("_begin called(already):".count($this->_dbs));
		$db->beginTransaction();
		$this->_dbs[]	= $db;
	}
	/**
	 * データベースコミット処理を行う．
	 * @param Zend_Db_Adapter_Mysqli $db
	 */
	protected function _commit() {
		$this->_log->debug("_commit called:".count($this->_dbs));
		foreach ($this->_dbs as $db){
			$db->commit();
		}
	}
	/**
	 * データベースロールバック処理を行う．
	 * @param Zend_Db_Adapter_Mysqli $db
	 */
	protected function _rollBack() {
		$this->_log->debug("_rollBack called:".count($this->_dbs));
		foreach ($this->_dbs as $db){
			$db->rollBack();
		}
	}
}