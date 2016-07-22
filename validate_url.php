<?php

$url = $_GET['url'];
$name = $_GET['nome'];
$link = $_GET['link'];
exec("python validator.py $url $name $link");

?>