#!/usr/bin/perl
use strict;
use MIME::Lite;
use Net::SMTP;
# This script will create a quick backup of the whole QMcBeaver directory
# simply call it, in a directory where there's a QMcBeaver directory with
# an optional message. e.g.
# ./makebackup.pl this is a message
#
# On cygwin, Net::SMTP was already set up
# I went here: http://search.cpan.org/dist/MIME-Lite/ for MIME::Lite
# Please change the email addresses... I don't want to recieve anybody else's backups! :-)

if(!(-d "QMcBeaver")){
  print "Run this in the same directory as a QMcBeaver directory!\n";
  die;
}

my $d = qx! date +%F.%H-%M-%S !;
chomp($d);
my $tarcommand = "tar --exclude CVS --exclude obj_* --exclude *.obj -cvzf";
my $basename = "QMcBeaver$d";
my $archivename = "$basename.tar.gz";
my $message = "@ARGV";
my $backupdir = "~/backups/";
print "Creating $archivename...\n";
print "$tarcommand $archivename QMcBeaver\n";
print "message: @ARGV\n";
`$tarcommand $archivename QMcBeaver`;

open (MESSENGER, ">$basename.txt");
print MESSENGER $message;
close(MESSENGER);

my $from = 'Amos Anderson <amosa@caltech.edu>';
my $to1 = 'nitroamos@gmail.com';
my $to2 = 'nitroamos@yahoo.com';
my $mail_host = 'smtp-server.its.caltech.edu';
my $subject = "Backing up $basename...";

print "Sending $archivename to $to1 and $to2...\n";
my $msg = MIME::Lite->new (
  From => $from,
  To => $to1,
  To => $to2,
  Subject => $subject,
  Type =>'multipart/mixed'
) or die "Error creating multipart container: $!\n";
$msg->attach (
  Type => 'TEXT',
  Data => $message,
) or die "Error adding the text message part: $!\n";
$msg->attach (
   Type => 'application/gzip',
   Path => $archivename,
   Filename => $archivename,
   Disposition => 'attachment'
) or die "Error adding $archivename: $!\n";
MIME::Lite->send('smtp', $mail_host, Timeout=>60);

#now the actual backing up part
$msg->send;
`mv $basename* $backupdir`;