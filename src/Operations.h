/*
 * Operations.h
 *
 *  Created on: Oct 26, 2016
 *      Author: mike
 */

#ifndef OPERATIONS_H_
#define OPERATIONS_H_

namespace ast {
struct Mode {};
struct Place : public Mode {};
struct AddIn : public Mode {};
struct SubIn : public Mode {};
struct RMult : public Mode {};
}

#endif /* OPERATIONS_H_ */
