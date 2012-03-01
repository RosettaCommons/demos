#!/usr/bin/perl
# return index of $1 in line of input (indexing from one)

#index of $param1 in param2, with first index == 1, or 0 if
sub indexArray_fromOne{
 1 while $_[0] ne pop;
 $#_+1
}

$line=<STDIN>;
@line_arr=split(/\s+/, $line);
print $line_arr[1], "\n";
print indexArray_fromOne(shift, @line_arr), "\n";
# print `'echo $line | sed ''s/\W\+/\n/g'' | cat -n | awk ''$2=="rmsALL"'' | awk ''{print $1}'' '`
