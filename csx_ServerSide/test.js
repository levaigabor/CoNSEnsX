var express = require('express');
var bodyParser = require('body-parser');
var fs = require("fs");
var util = require("util");
var app = express();

app.use(bodyParser({
        uploadDir: __dirname + 'uploads',
        keepExtensions: true
    }))

app.use(express.static(__dirname + '/public'));


// This route receives the posted form.
// As explained above, usage of 'body-parser' means
// that `req.body` will be filled in with the form elements
app.post('/', function(req, res, next){
  var userName = req.body.userName;
  var KARPLUS = req.body.KARPLUS;
  var superimpose = req.body.superimpose;
  var pdb_web = req.body.pdb_web;

  var pdb_upload = req.body.pdb_upload;






  var html = 'Hello: ' + userName + KARPLUS + superimpose +pdb_web+ '.<br>' +
             'PDB is: ' + pdb_upload +
             '<a href="/">Try again.</a>';
  res.send(html);

  console.log(req.body);
    console.log(req.files);
});

app.listen(8080);
