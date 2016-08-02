<?php

$url = $_GET['url'];
$name = $_GET['nome'];
$link = $_GET['link'];
exec("python tester.py $url $name $link");

?>