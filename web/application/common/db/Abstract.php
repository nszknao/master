<?php
abstract class
	DbAbstract
extends
	Zend_Db_Table_Abstract
{
	protected $_log	= null;
	
	/**
	 *
	 * ログ初期化.
	 */
	public function _loginit($className) {
		$log		= Zend_Registry::get('logger');
		$this->_log	= $log::getLogger($className);
		return	$this->_log;
	}
	/**
	 * 必要なカラム情報を含むパラメータのみ取り出す．
	 * @param $item
	 * @return array $result
	 */
	protected function setColumn($item) {
		$result	= array();
		$col	= $this->info();
		foreach ($col['cols'] as $colName) {
			if (array_key_exists($colName, $item)) {
				$result[$colName]	= $item[$colName];
			}
		}
		return $result;
	}
}