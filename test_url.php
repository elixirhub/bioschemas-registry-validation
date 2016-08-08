<?php

$url = $_GET['url'];
$name = $_GET['nome'];
$link = $_GET['link'];
$types = $_GET['types'];
exec("python tester.py $url $name $link $types");

?>