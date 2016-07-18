<?php
require_once 'Logger.php';

class Bootstrap extends Zend_Application_Bootstrap_Bootstrap
{
    /**
     *
     * Bootstrap Smarty view
     */
	protected function _initView() {
		$smartyConfig	= $this->getOption('smarty');
		$view = new Ext_View_Smarty($smartyConfig);
		$viewRenderer = Zend_Controller_Action_HelperBroker::getStaticHelper('ViewRenderer');
		$viewRenderer->setViewSuffix('tpl');
		$viewRenderer->setView($view);
		Zend_Registry::set('smarty', $smartyConfig);
		return $view;
	}
	/**
	 *
     * Bootstrap Logger
     * @return Logger
     */
	protected function _initLog() {
		$options = $this->getOptions();
		Logger::configure($options['log4php']['config']);
		define('LOGGER', 'default');
		$log	= Logger::getLogger('master-thesis');
		$log->debug("start");
		Zend_Registry::set('logger', $log);
		return $log;
	}
    /**
     *
     * Bootstrap db
     */
	protected function _initDb() {
		try {
			$options = new Zend_Config($this->getOptions('multidb'));
			$dbConfig = $options->multidb->db1;
			$dbConect = Zend_Db::factory($dbConfig->adapter, $dbConfig->params);
			Zend_Registry::set('db1', $dbConect);

			$dbConfig = $options->multidb->db2;
			$dbConect = Zend_Db::factory($dbConfig->adapter, $dbConfig->params);
			Zend_Registry::set('db2', $dbConect);
		} catch(Exception $e) {
			// データベース定義エラー
			syslog(LOG_ERR, "database Error:".$e->getMessage());
		}
		return;
	}
}