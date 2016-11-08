/*
 * Operations.h
 *
 *  Created on: Oct 26, 2016
 *      Author: mike
 */

#ifndef OPERATIONS_H_
#define OPERATIONS_H_

namespace ast {
/*
 * @class Place
 *
 * An empty type for emplacement assembly (replacing elements of a matrix
 * according to an assembler).
 */
struct Place {};

/*
 * @class AddIn
 *
 * An empty type for addition assembly (adding into elements of a matrix
 * according to an assembler)
 */
struct AddIn {};

/*
 * @class SubIn
 *
 * An empty type for subtraction assembly (subtracting from elements of a matrix
 * according to an assembler)
 */
struct SubIn {};

/*
 * @class RMult
 *
 * An empty type for performing right-multiplication according to an assembler
 */
struct RMult {};
}

#endif /* OPERATIONS_H_ */
