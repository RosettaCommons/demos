###
###
### This file is part of the CS-Rosetta Toolbox and is made available under
### GNU General Public License
### Copyright (C) 2011-2012 Oliver Lange
### email: oliver.lange@tum.de
### web: www.csrosetta.org
###
### This program is free software: you can redistribute it and/or modify
### it under the terms of the GNU General Public License as published by
### the Free Software Foundation, either version 3 of the License, or
### (at your option) any later version.
###
### This program is distributed in the hope that it will be useful,
### but WITHOUT ANY WARRANTY; without even the implied warranty of
### MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
### GNU General Public License for more details.
###
### You should have received a copy of the GNU General Public License
### along with this program.  If not, see <http://www.gnu.org/licenses/>.
###
###
#include <Python.h>

char *mmcif_get_string(void);

FILE *fp;

static PyObject *MMCIFlex_open_file(PyObject *self, PyObject *args)
{
	char *filename;

	if (!PyArg_ParseTuple(args, "s", &filename))
		return NULL;

	fp=fopen(filename, "r");	

	mmcif_set_file(fp);

	Py_INCREF(Py_None);

	return Py_None;
}	


static PyObject *MMCIFlex_close_file(PyObject *self, PyObject *args)
{
	/* verify no arguments */
	if (!PyArg_ParseTuple(args, ""))
		return NULL;

	fclose(fp);

	Py_INCREF(Py_None);

	return Py_None;
}	


static PyObject *MMCIFlex_get_token(PyObject *self, PyObject *args)
{
	int flag;
	char *value="";

	/* get token number */
	flag=mmcif_get_token();

	/* if flag==0 we are EOF */
	if(flag)
	{
		value=mmcif_get_string();
	}	

	/* return the (tokennumber, string) tuple */
	return Py_BuildValue("(is)", flag, value);
}


static PyMethodDef MMCIFlexMethods[]=
{
	{"open_file",	MMCIFlex_open_file, 	METH_VARARGS},
	{"close_file",	MMCIFlex_close_file,	METH_VARARGS},
	{"get_token",  	MMCIFlex_get_token, 	METH_VARARGS},
	{NULL,      	NULL}        			/* Sentinel */
};


void initMMCIFlex()
{
	(void) Py_InitModule("MMCIFlex", MMCIFlexMethods);
}
