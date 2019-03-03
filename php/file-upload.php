<?php
$timestamp = time();

$control_count = count($_FILES['control']['name']);
$admission_count = count($_FILES['admission']['name']);

exec("mkdir ../uploadFiles/" . $timestamp);
exec("mkdir ../uploadFiles/" . $timestamp . "/control");

for( $i=0 ; $i < $control_count ; $i++ ) {
  $tmpFilePath = $_FILES['control']['tmp_name'][$i];
  if ($tmpFilePath != ""){
    $newFilePath = "../uploadFiles/" . $timestamp . "/control/" . $_FILES['control']['name'][$i];
    if(move_uploaded_file($tmpFilePath, $newFilePath)) {
    }
  }
}
exec("gunzip ../uploadFiles/" . $timestamp . "/control/*.gz");

exec("mkdir ../uploadFiles/" . $timestamp . "/admission");
for( $i=0 ; $i < $admission_count ; $i++ ) {
  $tmpFilePath = $_FILES['admission']['tmp_name'][$i];
  if ($tmpFilePath != ""){
    $newFilePath = "../uploadFiles/" . $timestamp . "/admission/" . $_FILES['admission']['name'][$i];
    if(move_uploaded_file($tmpFilePath, $newFilePath)) {
    }
  }
}
exec("gunzip ../uploadFiles/" . $timestamp . "/admission/*.gz");
exec("find ../uploadFiles -type d -exec chmod 777 {} \;");
exec("screen -d -m -S ". $timestamp ." Rscript --vanilla ../r/chip_process.r ". dirname(getcwd())  ."/uploadFiles/" . $timestamp . " ". dirname(getcwd())  ."/uploadFiles/" . $timestamp . "/result.csv ". ($_POST["type"]=="affy" ? "affy" : "oligo")  ." > ../uploadFiles/" . $timestamp . "/outputFile.Rout 2>&1");
//exec("Rscript --vanilla ../r/chip_process.r ". dirname(getcwd())  ."/uploadFiles/" . $timestamp . " ". dirname(getcwd())  ."/uploadFiles/" . $timestamp . "/result.csv ". ($_POST["type"]=="affy" ? "affy" : "oligo")  ." > ../uploadFiles/" . $timestamp . "/outputFile.Rout 2>&1");
echo "<a href='/php/get-result.php?id=". $timestamp ."'>Link to your result</a>";?>
