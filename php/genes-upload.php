<html>
<head>
<style>
  img {
    float: left;
    width: 33%;
  }
</style>
</head>
<body>
<?php
$timestamp = time();
exec("mkdir ../uploadFiles/" . $timestamp);
$tmpFilePath = $_FILES['deg']['tmp_name'];
$newFilePath = "../uploadFiles/" . $timestamp . "/" . $_FILES['deg']['name'];
if ($tmpFilePath != ""){
  if(move_uploaded_file($tmpFilePath, $newFilePath)) {
  }
}
exec("python3 ../python/network_analysis.py \"". $_POST["gene_set"] ."\" \"". $newFilePath ."\"", $results);
print_r($results[0]);
echo "<br/>";
print_r($results[1]);
echo "<img src=\"../uploadFiles/" . $timestamp . "/compulsory.png\">";
echo "<img src=\"../uploadFiles/" . $timestamp . "/all.png\">";
echo "<img src=\"../uploadFiles/" . $timestamp . "/compulsory_additional.png\">";
?>
</body>
</html>
