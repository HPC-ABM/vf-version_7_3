/* 
 * File:  FieldVariable.h
 *
 * File Contents: Contains declarations for FieldVariable class
 *
 * Author: Yvonna
 * Contributors: Caroline Shung
 *               Nuttiiya Seekhao
 *               Kimberley Trickey
 *
 * Created on Jun 19, 2013, 9:58 PM
 *
 *****************************************************************************
 ***  Copyright (c) 2013 by A. Najafi-Yazdi                                ***
 *** This computer program is the property of Alireza Najafi-Yazd          ***
 *** and may contain confidential trade secrets.                           ***
 *** Use, examination, copying, transfer and disclosure to others,         ***
 *** in whole or in part, are prohibited except with the express prior     ***
 *** written consent of Alireza Najafi-Yazdi.                              ***
 *****************************************************************************/  // TODO(Kim): Update the file comment once we figure out the copyright issues

#ifndef FIELDVARIABLE_H
#define	FIELDVARIABLE_H

#include "../enums.h"

class World;

/*
 * FIELDVARIABLE CLASS DESCRIPTION:
 * FieldVariable is a parent class of WHChemical. 
 */
class FieldVariable {
 public:

    /*
     * Description:	Default FieldVariable constructor. 
     *
     * Return: void
     *
     * Parameters: void
     */
    FieldVariable();

    /*
     * Description:	FieldVariable constructor. 
     *
     * Return: void
     *
     * Parameters: orig  -- Pointer to an original FieldVariable
     */
    FieldVariable(const FieldVariable& orig);

    /*
     * Description:	Virtual FieldVariable destructor.
     *
     * Return: void
     *
     * Parameters: void
     */
    virtual ~FieldVariable();
};

#endif	/* FIELDVARIABLE_H */

