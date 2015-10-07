<!DOCTYPE html>
<html>
<head>
<meta charset="UTF-8">
<title>Title of the document</title>
</head>

<body>

<?php echo '<p>much PHP, such wow</p>'; ?>

<?php

    function crypto_rand_secure($min, $max)
    {
        $range = $max - $min;
        if ($range < 1) return $min; // not so random...
        $log = ceil(log($range, 2));
        $bytes = (int) ($log / 8) + 1; // length in bytes
        $bits = (int) $log + 1; // length in bits
        $filter = (int) (1 << $bits) - 1; // set all lower bits to 1
        do {
            $rnd = hexdec(bin2hex(openssl_random_pseudo_bytes($bytes)));
            $rnd = $rnd & $filter; // discard irrelevant bits
        } while ($rnd >= $range);
        return $min + $rnd;
    }

    function getToken($length)
    {
        $token = "";
        $codeAlphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
        $codeAlphabet.= "0123456789";
        $max = strlen($codeAlphabet) - 1;
        for ($i=0; $i < $length; $i++) {
            $token .= $codeAlphabet[crypto_rand_secure(0, $max)];
        }
        return $token;
    }

    function check_file_uploaded_name($filename)
    {
        (bool) ((preg_match("`^[-0-9A-Z_\.]+$`i",$filename)) ? true : false);
    }

    function check_file_uploaded_length($filename)
    {
        return (bool) ((mb_strlen($filename,"UTF-8") > 225) ? true : false);
    }

    $uploads_dir = '/var/www/consensx.itk.ppke.hu/public_html/file_uploads';
    $command = "consensx/consensx.py";

    // get ID for new calculation
    $csx_id = getToken(6);
    $command .= " -i " . $csx_id;
    echo "calculation id: " . $csx_id . "<br>";

    // upload PDB file if any
    if (isset($_FILES['pdb_upload'])) {
        $tmp_name = $_FILES["pdb_upload"]["tmp_name"];
        $PDB_name = $_FILES["pdb_upload"]["name"];

        $name_chars_ok   = check_file_uploaded_name($PDB_name);
        $name_size_is_ok = check_file_uploaded_length($PDB_name);

        echo "name is ". $name_chars_ok. "<br>";
        echo "size is ". $name_size_is_ok. "<br>";

        // if ($name_chars_ok and $name_size_is_ok) {
            move_uploaded_file($tmp_name, "$uploads_dir/$PDB_name");
            $command .= " -f " . "$uploads_dir/$PDB_name";
        // }
    }

    // upload STR file
    if (isset($_FILES['bmrb_upload'])) {
        $tmp_name = $_FILES["bmrb_upload"]["tmp_name"];
        $STR_name = $_FILES["bmrb_upload"]["name"];

        $name_chars_ok   = check_file_uploaded_name($STR_name);
        $name_size_is_ok = check_file_uploaded_length($STR_name);

        //if ($name_chars_ok and $name_size_is_ok) {
            move_uploaded_file($tmp_name, "$uploads_dir/$STR_name");
            $command .= " -b " . "$uploads_dir/$STR_name";
        //}
    }

    // upload NOE file if any
    if (isset($_FILES['xplor_upload'])) {
        echo "NOE FILE IS: " . $_FILES['xplor_upload'];
        $tmp_name = $_FILES["xplor_upload"]["tmp_name"];
        $NOE_name = $_FILES["xplor_upload"]["name"];

        $name_chars_ok   = check_file_uploaded_name($NOE_name);
        $name_size_is_ok = check_file_uploaded_length($NOE_name);

        // if ($name_chars_ok and $name_size_is_ok) {
            move_uploaded_file($tmp_name, "$uploads_dir/$NOE_name");
            $command .= " -r " . "$uploads_dir/$NOE_name";
        // }
    }

    // add fitting selection
    if (isset($_POST['superimpose'])) {
        echo "FITTING is enabled<br>";
        $command .= " -s";

        if (isset($_POST['fit_range'])) {
            echo 'FITTING range is ' . $_POST['fit_range'] . '<br>';
            $command .= " --fit_range " . $_POST['fit_range'];
        }
    }

    // add Karplus parameter set selection
    $command .= " -d ".$_POST['KARPLUS'];

    // add RDC SVD selection
    if (isset($_POST['RDCSVD'])) {
        echo "RDC is enabled<br>";
        $command .= " -R";
    }

    // add RDC LC model selection RDCLC
    $command .= " -l ". $_POST['RDCLC'];

    // add r^-3 averaging selection
    if (isset($_POST['r3average'])) {
        echo "r<sup>3</sup> averaging is enabled<br>";
        $command .= " -r3";
    }

    echo $command;
    // $output = shell_exec('python3 consensx/consensx.py -b consensx/dummy.str -f consensx/1d3z.pdb -i ' . $csx_id);
    $output = shell_exec('python3 ' . $command);

    echo "<br>" . $output . "<br>";

    // header("Location: http://consensx.itk.ppke.hu/calculations/" . $csx_id . "/result_sheet.html");
    // exit;


    // $dir = 'consensx/myDir';

     // create new directory with 744 permissions if it does not exist yet
     // owner will be the user/group the PHP script is run under
     // if ( !file_exists($dir) ) {
     //  mkdir ($dir, 0764);
     // }

     // file_put_contents ($dir.'/test.txt', 'Hello File');
?>

</body>
</html>

