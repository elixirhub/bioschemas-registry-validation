<?php

$url = $_GET['url'];
$name = $_GET['nome'];
$link = $_GET['link'];
$types = $_GET['types'];
$mail = $_GET['mail'];
exec("python validator.py $url $name $link $types $mail");

?>