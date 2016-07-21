<?php

$url = $_GET['url'];
$name = $_GET['nome'];
exec("python validator.py $url $name");

?>