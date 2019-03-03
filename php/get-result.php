<?php
if (isset($_GET["id"]) && file_exists("../uploadFiles/" . $_GET["id"] . "/result.csv")) {
  $f = fopen("../uploadFiles/" . $_GET["id"] . "/result.csv", "r");?>
  <html>
    <head>
      <title>Result table</title>
      <meta charset="utf-8">
      <meta name="viewport" content="width=device-width, initial-scale=1">
      <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.4.0/css/bootstrap.min.css">
      <script src="https://ajax.googleapis.com/ajax/libs/jquery/3.3.1/jquery.min.js"></script>
      <script src="https://maxcdn.bootstrapcdn.com/bootstrap/3.4.0/js/bootstrap.min.js"></script>
      <script src="https://cdn.datatables.net/1.10.19/js/jquery.dataTables.min.js"></script>
      <script type="text/javascript">
        $(document).ready(function () {
          $('#dtBasicExample').DataTable();
          $('.dataTables_length').addClass('bs-select');
        });
      </script>
      <style>
        table.dataTable thead .sorting:after,
        table.dataTable thead .sorting:before,
        table.dataTable thead .sorting_asc:after,
        table.dataTable thead .sorting_asc:before,
        table.dataTable thead .sorting_asc_disabled:after,
        table.dataTable thead .sorting_asc_disabled:before,
        table.dataTable thead .sorting_desc:after,
        table.dataTable thead .sorting_desc:before,
        table.dataTable thead .sorting_desc_disabled:after,
        table.dataTable thead .sorting_desc_disabled:before {
          bottom: .5em;
        }
      </style>
      </head>
      <body>
	<?php echo "Click to download <a href=\"../uploadFiles/". $_GET["id"] . "/result.csv\">full results csv</a> or <a href=\"../uploadFiles/". $_GET["id"] . "/result_filtered.csv\">filtered results csv</a>"?>
        <table id="dtBasicExample" class="table table-striped table-bordered table-sm" cellspacing="0" width="100%">
          <thead>
          <?php
            $line = fgetcsv($f);
            echo "<tr>";
            foreach ($line as $cell) {
              echo "<th>" . htmlspecialchars($cell) . "</th>";
            }
            echo "</tr>\n";
            ?>
          </thead>
        <tbody>
        <?php
          while (($line = fgetcsv($f)) !== false) {
            echo "<tr>";
            foreach ($line as $cell) {
              $value = htmlspecialchars($cell);
              echo "<td>" . sprintf((is_numeric($value) ? (strval($value) !== strval(intval($value)) ? "%f" : "%d") : "%s"), $value) . "</td>";
              //echo "<td>" . sprintf("%f", htmlspecialchars($cell)) . "</td>";
            }
            echo "</tr>\n";
          }
          ?>
        </tbody>
      </table>
    </body>
  </html>
  <?php
  fclose($f);
} else {
  echo "Please wait, your data is being processed<br/>";
  echo "<a href='/php/get-result.php?id=". $_GET["id"] ."'>Reload link</a>";
}
?>
