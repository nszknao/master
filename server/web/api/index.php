<?php
function myloader($classname) {
	return @include str_replace('_', DIRECTORY_SEPARATOR, $classname).'.php';
}
spl_autoload_register('myloader');

// Define path to application directory
defined('APPLICATION_PATH')		|| define('APPLICATION_PATH', realpath(dirname(__FILE__) . '/../application'));
defined('APPLICATION_ENV')		|| define('APPLICATION_ENV', (getenv('APPLICATION_ENV') ? getenv('APPLICATION_ENV') : 'production'));
defined('APPLICATION_TYPE')		|| define('APPLICATION_TYPE', (getenv('APPLICATION_TYPE') ? getenv('APPLICATION_TYPE') : 'api'));

// Ensure library/ is on include_path
set_include_path(implode(PATH_SEPARATOR, array(
    realpath(APPLICATION_PATH . '/../library'),
    get_include_path(),
)));

// PHPの環境設定
require_once 'env/phpenv.php';

/** Zend_Application */
require_once 'Zend/Application.php';

// Smarty
require_once 'Smarty/Smarty.class.php';

// Create application, bootstrap, and run
$application = new Zend_Application(
    APPLICATION_ENV,
    APPLICATION_PATH . '/../configs/application.ini'
);
$application->bootstrap()
            ->run();