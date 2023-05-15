#include <iostream>
#include <vector>
#include <cmath>
using namespace std;

// ################## REGIONES ################################
vector< vector<float> > allSizes;

// ################## CLASE PUNTO ##############################
// DIMENSIONES DEL POKEMON
template <int ndimension>
class Point{
    private:
        float datosPokemon[ndimension]; // height_m, hp, weight_kg_std

    public:
        // ################# CONSTRUCTOR ###########################
        Point(){
            for(int i = 0; i < ndimension; i++){
                datosPokemon[i] = 0;
            }
        }

        Point( vector<float> _datosPokemon){
            for(int i = 0; i < ndimension; i++){
                datosPokemon[i] = _datosPokemon[i];
            }
        }

        // ################# METODOS ###########################
        float  &operator[](int posicion){
            return datosPokemon[posicion];
        }
};

// ################## CLASE POKEMON ##############################
// COORDENADAS DEL POKEMON
template <int ndimension>
class Pokemon{
    private:
        float cordsPokemon[ndimension];

    public:
        // ################# CONSTRUCTOR ###########################
        Pokemon(){
            for(int i = 0; i < ndimension; i++){
                cordsPokemon[i] = 0;
            }
        }

        Pokemon( vector<float> _cordsPokemon ){
            for(int i = 0; i < ndimension; i++){
                cordsPokemon[i] = _cordsPokemon[i];
            }
        }

        // ################# METODOS ###########################
        float  &operator[](int posicion){
            return cordsPokemon[posicion];
        }
};

// ################## CLASE RECORD ##############################
// ALL DATA OF POKEMON
template< typename T, int ndim >
class Record {
    private:
        T object; // Coordenadas del pokemon
        Point<ndim> tuple; //  Dimensiones del pokemon

    public:
        // ################# CONSTRUCTOR ###########################
        Record(float attack = 0, float defense = 0, float speed = 0, float height_m = 0, float hp = 0, float weight_kg_std = 0){

            // COORDENADAS
            vector<float> auxCoordenadas;
            auxCoordenadas.push_back(attack);
            auxCoordenadas.push_back(defense);
            auxCoordenadas.push_back(speed);
            Pokemon<ndim> _object(auxCoordenadas);
            object = _object;

            // DIMENSIONES
            vector<float> auxDimensiones;
            auxDimensiones.push_back(height_m);
            auxDimensiones.push_back(hp);
            auxDimensiones.push_back(weight_kg_std);
            Point<ndim> _tuple(auxDimensiones);
            tuple = _tuple;
        }

        // ################# METODOS ###########################
        T getObject(){
            return  object;
        }

        Point<ndim> getTupla(){
            return tuple;
        };
};

// ################## CLASE NODO ##############################
template< typename T, int ndim >
class Node {
    private:
        vector< Record<T,ndim> > grupoPokemones;    // Vector of records
        vector< Node<T,ndim> > grupoHijos;          // Vector of childrens type Node
        Point<ndim> bottomLeft;                     // Point x,y,z
        Point<ndim> upperRight;                     // Point x,y,z
        int maxRecords;                             // Limit of records per node
        int countRecords;                           // Data counter per node
        bool isLeaf;                                // Flag

    public:
        // ################# CONSTRUCTOR ###########################
        Node(){
            maxRecords = 0;
            countRecords = 0;
            isLeaf = true;

            bottomLeft = Point<ndim>();
            upperRight = Point<ndim>();
        }

        Node(int _maxhijos){
            maxRecords = _maxhijos;
            countRecords = 0;
            isLeaf = true;

            bottomLeft = Point<ndim>();
            upperRight = Point<ndim>();
        }

        // ################# METODOS ###########################

        // ----------- Functions to print ------------
        void printRecord(Node *pntN){
            cout << "Numero de records guardados: " << pntN->grupoPokemones.size() << endl;
            for(int j = 0; j < pntN->grupoPokemones.size(); j++){
                Pokemon<ndim> getObject = pntN->grupoPokemones[j].getObject();
                Point<ndim> getTupla = pntN->grupoPokemones[j].getTupla();

                cout << "[ ";
                for(int i = 0; i < ndim; i++){
                    cout << getObject[i] << " - ";
                }
                for(int i = 0; i < ndim; i++){
                    cout << getTupla[i] << " - ";
                }
                cout << " ]" << endl;
            }
        }

        void printRecord2(Record<T, ndim> record){
            Pokemon<ndim> getObject = record.getObject();
            Point<ndim> getTupla = record.getTupla();
            cout << "[ ";
            for(int i = 0; i < ndim; i++){
                cout << getObject[i] << " - ";
            }
            for(int i = 0; i < ndim; i++){
                cout << getTupla[i] << " - ";
            }
            cout << " ]" << endl;
        }

        void printBottomLeft(Node *pntN){
            cout << "BottomLeft Actual es: " << endl;
            for(int i = 0; i < ndim; i++){
                cout << pntN->bottomLeft[i] << " - ";
            }
        }

        void printUpperRight( Node *pntN ){
            cout << endl << "UpperRight Actual es: " << endl;
            for(int i = 0; i < ndim; i++){
                cout << pntN->upperRight[i] << " - ";
            }
        }

        // ----------- Functions to recalculate and get MBR ------------
        Point<ndim> recalculateBottomLeft(Node *pntN){
            Point<ndim> result = Point<ndim>();

            for(int j = 0; j < ndim; j++){
                Pokemon<ndim> getObjectAux1 = pntN->grupoPokemones[0].getObject();
                float minNum1 = getObjectAux1[j];
                for(int i = 1; i < pntN->countRecords; i++){
                    Pokemon<ndim> getObjectAux2 = pntN->grupoPokemones[i].getObject();
                    float minNum2 = getObjectAux2[j];
                    if(minNum2 < minNum1){
                        minNum1 = minNum2;
                    }
                }
                result[j] = minNum1;
            }
            return result;
        }

        Point<ndim> recalculateUpperRight(Node *pntN){
            //Pokemon<ndim> getObject = pntN->grupoPokemones[0].getObject();
            //Point<ndim> getTupla = pntN->grupoPokemones[0].getTupla();
            /*
            cout << "Cords :";
            for(int i = 0; i < 3; i++){
                cout << getObject[i] << " - ";
            }
            cout << endl << "Dimensiones :";
            for(int i = 0; i < 3; i++){
                cout << getTupla[i] << " - ";
            }
            cout << endl;
            */

            Point<ndim> result = Point<ndim>();

            for(int j = 0; j < ndim; j++){
                Point<ndim> getTuplaAux1 = pntN->grupoPokemones[0].getTupla();
                Pokemon<ndim> getObjectAux1 = pntN->grupoPokemones[0].getObject();
                float maxNum1 =  getObjectAux1[j] + getTuplaAux1[j];
                for(int i = 1; i < pntN->countRecords; i++){
                    Point<ndim> getTuplaAux2 = pntN->grupoPokemones[i].getTupla();
                    Pokemon<ndim> getObjectAux2 = pntN->grupoPokemones[i].getObject();
                    float maxNum2 = getObjectAux2[j] + getTuplaAux2[j];
                    if(maxNum2 > maxNum1){
                        maxNum1 = maxNum2;
                    }
                }
                result[j] = maxNum1;
            }
            return result;
        }

        pair< Point<ndim>, Point<ndim> > get_BL_UR_Record( Record<T, ndim> record ){

            pair< Point<ndim>, Point<ndim> > result;

            // Sacamos el BottomLeft del nuevo record a insertar
            Pokemon<ndim> getObjectAux1 = record.getObject();
            vector<float> auxBLRecord;
            for(int z = 0; z < ndim; z++){
                auxBLRecord.push_back( getObjectAux1[z] );
            }
            Point<ndim> BLRecord(auxBLRecord);
            result.first = BLRecord;

            // Sacamos el UpperRight del nuevo record a insertar
            Point<ndim> getTuplaAux1 = record.getTupla();
            vector<float> auxURRecord;
            for(int z = 0; z < ndim; z++){
                auxURRecord.push_back( getObjectAux1[z] + getTuplaAux1[z] );
            }
            Point<ndim> URRecord(auxURRecord);
            result.second = URRecord;

            return result;
        }

        // ----------- Functions to PeekSeeds ------------
        float calculateArea(Point<ndim> bottomLeft, Point<ndim> upperRight){
            float area = 1;
            for(int i = 0; i < ndim; i++){
                float aux = upperRight[i] - bottomLeft[i];
                area = area * aux;
            }

            return area;
        }

        pair< Point<ndim>, Point<ndim> > getDimensionsJ( Record<T, ndim>  E1I, Record<T, ndim> E2I){

            pair<Point<ndim>, Point<ndim>> BL_UR_Record1 = get_BL_UR_Record(E1I);
            Point<ndim> BLRecord1 = BL_UR_Record1.first;
            Point<ndim> UPRecord1 = BL_UR_Record1.second;

            pair<Point<ndim>, Point<ndim>> BL_UR_Record2 = get_BL_UR_Record(E2I);
            Point<ndim> BLRecord2 = BL_UR_Record2.first;
            Point<ndim> UPRecord2 = BL_UR_Record2.second;

            pair< Point<ndim>, Point<ndim> > result;

            vector<float> resAuxBL;
            for(int i = 0; i < ndim; i++){
                float aux1 = BLRecord1[i];
                float aux2 = BLRecord2[i];
                if(aux2 < aux1){
                    aux1 = aux2;
                }
                resAuxBL.push_back(aux1);
            }
            Point<ndim> auxResultBL(resAuxBL);

            vector<float> resAuxUR;
            for(int i = 0; i < ndim; i++){
                float aux1 = UPRecord1[i];
                float aux2 = UPRecord2[i];
                if(aux2 > aux1){
                    aux1 = aux2;
                }
                resAuxUR.push_back(aux1);
            }
            Point<ndim> auxResultUR(resAuxUR);

            result.first = auxResultBL;
            result.second = auxResultUR;

            /*
            cout << "----- TEST ------" << endl;
            cout << "Bottom Left of J: ";
            for(int i = 0; i < ndim; i++){
                cout << auxResultBL[i] << " - ";
            }
            cout << endl;

            cout << "Upper Right of J: ";
            for(int i = 0; i < ndim; i++){
                cout << auxResultUR[i] << " - ";
            }
            cout << endl;
            */

            return  result;
        }

        pair< Point<ndim>, Point<ndim> > getDimensionsE( Record<T, ndim>  E_NI){

            pair<Point<ndim>, Point<ndim>> BL_UR_Record1 = get_BL_UR_Record(E_NI);
            Point<ndim> BLRecord1 = BL_UR_Record1.first;
            Point<ndim> UPRecord1 = BL_UR_Record1.second;

            pair< Point<ndim>, Point<ndim> > result;

            result.first = BLRecord1;
            result.second = UPRecord1;

            /*
            cout << "----- TEST ------" << endl;
            cout << "Bottom Left of E: ";
            for(int i = 0; i < ndim; i++){
                cout << BLRecord1[i] << " - ";
            }
            cout << endl;

            cout << "Upper Right of E: ";
            for(int i = 0; i < ndim; i++){
                cout << UPRecord1[i] << " - ";
            }
            cout << endl;
            */

            return  result;
        }

        // ----------- Functions to PeekNext ------------
        float calculateEnlargment( Record<T, ndim> &record, Node *hijo){
            //  --> Sacamos el BottomLeft y UpperRight de la region/hijo[i]
            Point<ndim> BLRegion = hijo->bottomLeft;
            Point<ndim> URRegion = hijo->upperRight;

            // --> Sacamos el BottomLeft y UpperRight del nuevo record a insertar
            Pokemon<ndim> getObjectAux1 = record.getObject();
            vector<float> auxBLRecord;
            for(int z = 0; z < ndim; z++){
                auxBLRecord.push_back( getObjectAux1[z] );
            }
            Point<ndim> BLRecord(auxBLRecord);

            Point<ndim> getTuplaAux1 = record.getTupla();
            vector<float> auxURRecord;
            for(int z = 0; z < ndim; z++){
                auxURRecord.push_back( getObjectAux1[z] + getTuplaAux1[z] );
            }
            Point<ndim> URRecord(auxURRecord);

            // -------------- Print ----------------------
            /*
            cout << endl << "Region evaluando..." << endl;
            printBottomLeft(hijo);
            printUpperRight(hijo);
            cout << endl;

            cout << "Nuevo Record " << endl;
            cout << "BT: " << endl;
            for(int a = 0; a < ndim; a++){
                cout << BLRecord[a] << " - ";
            }
            cout << endl << "UR: " << endl;
            for(int a = 0; a < ndim; a++){
                cout << URRecord[a] << " - ";
            }
            cout << endl << endl;
             */
            // ------------------------------------------

            // --> Sacamos el area de la region actual
            float areaRegion = 1;
            for(int i = 0; i < ndim; i++){ // El area se la multiplicacion de todas las dimenciones
                float areaRegionAux = URRegion[i] - BLRegion[i];
                areaRegion = areaRegion * areaRegionAux;
            }

            // --> Sacamos el area de MBR expandido (con el nuevo record)
            float areaTotal = 1;
            for(int i = 0; i < ndim; i++){ // El area se la multiplicacion de todas las dimenciones

                float LimInf = min(BLRegion[i], BLRecord[i]);
                float LimSup = max(URRegion[i], URRecord[i]);

                float areaTotalAux = LimSup - LimInf;

                areaTotal = areaTotal * areaTotalAux;
            }

            // --> Restamos ambas para sacar el area de expancion
            float areaExpacion = areaTotal - areaRegion;
            //cout << "AreaTotal (" << areaTotal << ") - AreaRegion(" << areaRegion << ") = areaExpacion (" << areaExpacion << ")" << endl;
            return areaExpacion;
        }

        // ----------- Functions for Cuadratic Split ------------
        pair< Node, Node > peekSeeds_nextPeek( Node *pntN ){

            // ########################### PEEK SEEDS ################################################
            //cout << "Buscando semillas.... " << endl;
            Node<T, ndim> L; // root
            Node<T, ndim> *pntNL; // puntero
            Node<T, ndim> auxL(maxRecords);
            L = auxL;
            pntNL = &L;

            Node<T, ndim> LL; // root
            Node<T, ndim> *pntNLL; // Puntero
            Node<T, ndim> auxLL(maxRecords);
            LL = auxLL;
            pntNLL = &LL;

            // -------------- LOGICA ---------------------
            int numeroSemilla1 = 0;
            int numeroSemilla2 = 1;

            //cout << "Evaluando con punto " << numeroSemilla1 << " y " << numeroSemilla2 << endl;

            pair< Point<ndim>, Point<ndim> > areaJInicialAux = getDimensionsJ(pntN->grupoPokemones[numeroSemilla1], pntN->grupoPokemones[numeroSemilla2]);
            float areaJInicial = calculateArea(areaJInicialAux.first, areaJInicialAux.second);

            pair< Point<ndim>, Point<ndim> > areaE1InicialAux = getDimensionsE(pntN->grupoPokemones[numeroSemilla1]);
            float areaE1Inicial = calculateArea(areaE1InicialAux.first, areaE1InicialAux.second);

            pair< Point<ndim>, Point<ndim> > areaE2InicialAux = getDimensionsE(pntN->grupoPokemones[numeroSemilla2]);
            float areaE2Inicial = calculateArea(areaE2InicialAux.first, areaE2InicialAux.second);

            float dInicial = areaJInicial - areaE1Inicial - areaE2Inicial;
            //cout << " ---> Desperdicio acumulado: " << dInicial << endl;

            int j = 2;
            for(int i = 0; i < pntN->countRecords; i++){
                for( j ; j < pntN->countRecords; j++){
                    //cout << endl << "Evaluando con punto " << i << " con " << j << endl;

                    pair< Point<ndim>, Point<ndim> > areaJAux = getDimensionsJ(pntN->grupoPokemones[i], pntN->grupoPokemones[j]);
                    float areaJ = calculateArea(areaJAux.first, areaJAux.second);

                    pair< Point<ndim>, Point<ndim> > areaE1Aux = getDimensionsE(pntN->grupoPokemones[i]);
                    float areaE1 = calculateArea(areaE1Aux.first, areaE1Aux.second);

                    pair< Point<ndim>, Point<ndim> > areaE2Aux = getDimensionsE(pntN->grupoPokemones[j]);
                    float areaE2 = calculateArea(areaE2Aux.first, areaE2Aux.second);

                    float d = areaJ - areaE1 - areaE2;
                    //cout << " ---> Desperdicio acumulado: " << d << endl;

                    if( d > dInicial ){
                        numeroSemilla1 = i;
                        numeroSemilla2 = j;
                    }
                }
                j = i+2;
            }
            //cout << endl << "Los records mas lejanos son: Record " << numeroSemilla1 << " y Record " << numeroSemilla2 << endl;
            //cout << "Insertando semillas ..." << endl;

            insertNode( pntN->grupoPokemones[numeroSemilla1], pntNL );
            insertNode( pntN->grupoPokemones[numeroSemilla2], pntNLL );

            // ########################### NEXT PEEK ################################################
            //cout << "Insertando records faltantes .... " << endl;
            for(int i = 0; i < pntN->countRecords; i++){
                if( i != numeroSemilla1 && i != numeroSemilla2 ){
                    float crecimientoSeed1 = calculateEnlargment(pntN->grupoPokemones[i], pntNL);
                    float crecimientoSeed2 = calculateEnlargment(pntN->grupoPokemones[i], pntNLL);

                    if(crecimientoSeed1 < crecimientoSeed2){
                        //cout << endl << "Insertando en la semilla 1" << endl;
                        insertNode(pntN->grupoPokemones[i],pntNL);
                    }
                    else{
                        //cout << endl << "Insertando en la semilla 2" << endl;
                        insertNode(pntN->grupoPokemones[i], pntNLL);
                    }
                }
                else{
                    //cout << "El hijo " << i << " es una semilla y ya fue insertado" << endl;
                }
            }

            pair <Node, Node > resultPeekSeeds_NextPeek;
            resultPeekSeeds_NextPeek.first = L;
            resultPeekSeeds_NextPeek.second = LL;

            return resultPeekSeeds_NextPeek;
        }

        pair< Node, Node > splitCuadratico(Node *pntN){

            pair <Node , Node> resultSplit;

            Node<T, ndim> L; // root
            Node<T, ndim> *pntNL; // puntero
            Node<T, ndim> auxL(maxRecords);
            L = auxL;
            pntNL = &L;

            Node<T, ndim> LL; // root
            Node<T, ndim> *pntNLL; // Puntero
            Node<T, ndim> auxLL(maxRecords);
            LL = auxLL;
            pntNLL = &LL;

            // ---------------- PEEKSEEDS y NEXTPEEK --------------------------
            pair <Node, Node > resultPeekSeeds = peekSeeds_nextPeek(pntN);
            L = resultPeekSeeds.first;
            LL = resultPeekSeeds.second;

            //cout << " -------- RESULTADO TRAS SPLIT ---------" << endl;
            //printRecord(pntNL);
            //cout << endl;
            //printRecord(pntNLL);

            resultSplit.first = L;
            resultSplit.second = LL;

            return  resultSplit;
        }

        // ----------- Functions for Equidistant Split ------------
        float distance2Pokemons( Pokemon<ndim> pntBottomLeft, Point<ndim> pntUpperRigth ){
            float res = 0;
            for(int i = 0; i < ndim; i++){
                float aux = pntUpperRigth[i] - pntBottomLeft[i];
                res = res + pow(aux, 2);
            }
            return sqrt(res);
        }

        pair< Node, Node > splitEquisdistantes(Node *pntN){
            // ################ PEEK SEEDS #########################
            pair <Node , Node> resultSplit;

            //1. Seteamos los primeros 2 puntos (semillas)
            vector<float> smlla1, smlla2;
            for (int i = 0; i < ndim; i++){
                float mitad = (pntN->bottomLeft[i] + pntN->upperRight[i])/2 ;

                float equisAux1;
                equisAux1 = ( pntN->bottomLeft[i] + mitad )/2;
                smlla1.push_back(equisAux1);

                float equisAux2;
                equisAux2 = ( mitad + pntN->upperRight[i] )/2;
                smlla2.push_back(equisAux2);
            }

            cout << endl;

            /*
            cout << " --> Equisdistante 1: ";
            for(int i = 0; i < ndim; i++){
                cout << smlla1[i] << " - ";
            }
            cout << endl;
            cout << " --> Equisdistante 2: ";
            for(int i = 0; i < ndim; i++){
                cout << smlla2[i] << " - ";
            }
            cout << endl;
             */

            Point<ndim> equisdistante1(smlla1); // Semilla 1
            Point<ndim> equisdistante2(smlla2); // Semilla 2

            // ################# RE INSERT ########################
            Node<T, ndim> L; // root
            Node<T, ndim> *pntNL;
            Node<T, ndim> auxL(maxRecords);
            L = auxL;
            pntNL = &L;

            Node<T, ndim> LL; // root
            Node<T, ndim> *pntNLL;
            Node<T, ndim> auxLL(maxRecords);
            LL = auxLL;
            pntNLL = &LL;

            // Comparamos cada hijo del puntero con los equisdistantes que obtuvimos
            for(int i = 0; i < pntN->countRecords; i++){

                // Comparo hijo con los 2 equisdistantes
                float auxDistance1 = distance2Pokemons(pntN->grupoPokemones[i].getObject(), equisdistante1);
                float auxDistance2 = distance2Pokemons(pntN->grupoPokemones[i].getObject(), equisdistante2);

                cout << endl << "..... Insertando datos guardados " << endl;
                if(auxDistance1 < auxDistance2){
                    cout << "P" << i+1 << " va al equisdistante 1" << endl;
                    // Guardo en un nodo el dato insertado
                    cout << "Insertando" << endl;
                    insertNode(pntN->grupoPokemones[i], pntNL);
                    //auxNode1.insertNode(pntN->grupoPokemones[i], auxNode1);
                }
                else{
                    cout << "P" << i+1 << " va al equisdistante 2" << endl;
                    insertNode(pntN->grupoPokemones[i], pntNLL);
                    //auxNode2.insertNode(pntN->grupoPokemones[i], auxNode2);
                }
            }

            cout << " ---- RESULTADO TRAS SPLIT ---------" << endl;
            resultSplit.first = L;
            resultSplit.second = LL;

            printRecord(pntNL);
            cout << endl;
            printRecord(pntNLL);

            return  resultSplit;
        }

        // ----------- INSERT FUNCTION ----------------------
        bool insertNode(Record<T,ndim> &record, Node *pntN){
            //cout <<  endl << "INSERTANDO: ";
            //printRecord2(record);

            //cout << "-------- Preview ----------" << endl;
            //printRecord(pntN);
            //cout << "--------------------------" << endl;

            if(pntN->isLeaf == true) {
                if(pntN->countRecords < pntN->maxRecords) {
                    // --------------- INSERT NORMAL -----------------------
                    pntN->grupoPokemones.push_back(record);
                    pntN->countRecords = pntN->countRecords + 1;

                    // -->  Recalcular el MBR
                    pntN->bottomLeft = recalculateBottomLeft(pntN);
                    pntN->upperRight = recalculateUpperRight(pntN);

                    // --> Imprimir MBR
                    //printBottomLeft(pntN);
                    //printUpperRight(pntN);

                    //cout << endl << "------------- Se inserto correctamente ------------------" << endl;
                    return true;
                }

                else{
                    //cout << "Supera granularidad... -> Split()" << endl;
                    // ------------------ INSERT CON SPLIT -----------------------
                    pntN->countRecords = pntN->countRecords + 1;
                    pntN->grupoPokemones.push_back(record);

                    // -->  Recalcular el MBR
                    pntN->bottomLeft = recalculateBottomLeft(pntN);
                    pntN->upperRight = recalculateUpperRight(pntN);

                    // --> Imprimir MBR
                    //printBottomLeft(pntN);
                    //printUpperRight(pntN);
                    //cout << endl << endl;

                    // --> split()
                    pntN->isLeaf = false;
                    //pair<Node, Node> auxSplit = splitEquisdistantes(pntN);     // Usando equisdistantes
                    pair<Node, Node> auxSplit = splitCuadratico(pntN);           // Usando Metodo cuadratico

                    // --> Añadimos los hijos al nodo actual
                    pntN->grupoHijos.push_back(auxSplit.first);
                    pntN->grupoHijos.push_back(auxSplit.second);

                    return  true;
                }
            }
            else{
                //cout <<  "El root NO es un una hoja... RecorrerHijos() " << endl;
                // ------------------ INSERT A UN NODO HIJO ---------------------
                pntN->countRecords = pntN->countRecords + 1;
                pntN->grupoPokemones.push_back(record);

                // -->  Recalcular el MBR
                pntN->bottomLeft = recalculateBottomLeft(pntN);
                pntN->upperRight = recalculateUpperRight(pntN);

                // --> Imprimir MBR
                //printBottomLeft(pntN);
                //printUpperRight(pntN);

                // a. Verificar si esta dentro de un hijo
                bool isInMBR;
                for(int i = 0; i < pntN->grupoHijos.size(); i++){

                    isInMBR = true;

                    //  --> Obtenemos BottomLeft y UpperRight de cada hijo
                    Node<T, ndim> *pntNAux;
                    pntNAux = &pntN->grupoHijos[i];
                    Point<ndim> BLRegion = pntNAux->bottomLeft;
                    Point<ndim> URRegion = pntNAux->upperRight;

                    // --> Obtenemos BottomLeft y UpperRight del record a insertar
                    Pokemon<ndim> getObjectAux1 = record.getObject();
                    vector<float> auxBLRecord;
                    for(int z = 0; z < ndim; z++){
                        auxBLRecord.push_back( getObjectAux1[z] );
                    }
                    Point<ndim> BLRecord(auxBLRecord); // BLRecord = BottomLeft del nuevo record

                    Point<ndim> getTuplaAux1 = record.getTupla();
                    vector<float> auxURRecord;
                    for(int z = 0; z < ndim; z++){
                        auxURRecord.push_back( getObjectAux1[z] + getTuplaAux1[z] );
                    }
                    Point<ndim> URRecord(auxURRecord); // URRecord = UpperRight del nuevo record

                    // --> Imprimirmos hijo
                    /*
                    cout << endl << "Dimensiones de los hijos: " << endl;
                    printBottomLeft(pntNAux);
                    printUpperRight(pntNAux);
                    cout << endl;
                    */

                    // --> Imprimirmos record
                    /*
                    cout << "Dimensiones del nuevo record: " << endl;
                    cout << "BT: " << endl;
                    for(int a = 0; a < ndim; a++){
                        cout << BLRecord[a] << " - ";
                    }
                    cout << endl << "UR: " << endl;
                    for(int a = 0; a < ndim; a++){
                        cout << URRecord[a] << " - ";
                    }
                    cout << endl << endl;
                    */

                    // --> Operaciones para comprobar si esta dentro
                    int auxContador = 0;
                    for(int ejeN = 0; ejeN < ndim; ejeN++){
                        float limit1 = BLRegion[ejeN];
                        float limit2 = URRegion[ejeN];

                        //cout << "-> EJE " << ejeN << endl;
                        //cout << "Comparando: " << limit1 << " <= " << BLRecord[auxContador] << " <= " << limit2 << " y " << limit1 << " <= " << URRecord[auxContador] << " <= " << limit2 << endl;
                        if( limit1 <= BLRecord[auxContador] && BLRecord[auxContador] <= limit2 && limit1 <= URRecord[auxContador] && URRecord[auxContador] <= limit2 ){
                            //cout << "Continuamos..." << endl;
                        }
                        else{
                            //cout << "No esta dentro del hijo " << i+1 << ". Pasando al siguiente hijo ..."<< endl;
                            isInMBR = false;
                            break;
                        }
                        auxContador++;
                    }

                    if(isInMBR == true){
                        //cout << "SI esta dentro de un hijo ... Insertamos" << endl;
                        insertNode(record, pntNAux);
                        break;
                    }
                }

                // b. No esta dentro de un hijo. Ejecutar ChoosseLeaf
                if(isInMBR == false){
                    //cout << endl << "NO esta dentro de un hijo ... Ejecutar ChoosseLeaf"<< endl;

                    int posicionFinal = 0;

                    // --> Calcular crecimiento/expancion del primer hijo con record
                    Node<T, ndim> *pntChosseLeaf1;
                    pntChosseLeaf1 = &pntN->grupoHijos[0];
                    float crecimientoInicial = calculateEnlargment(record, pntChosseLeaf1);

                    // --> Comparamos con el crecimiento/expancion de todos los otros hijos
                    for(int i = 1; i < pntN->grupoHijos.size(); i++ ){
                        Node<T, ndim> *pntChosseLeaf2;
                        pntChosseLeaf2 = &pntN->grupoHijos[i];
                        float crecimientoAux = calculateEnlargment(record, pntChosseLeaf2);

                        // --> Guardamos la posicion del menor crecimiento en la variable "Posicion Final"
                        if(crecimientoAux < crecimientoInicial){
                            posicionFinal = i;
                            crecimientoInicial = crecimientoAux;
                        }
                    }

                    //cout << endl << "El hijo con la menor expancion es el numero " << posicionFinal << endl;

                    // --> Insertamos en el hijo correspondiente
                    Node<T, ndim> *pntChooseLeaf;
                    pntChooseLeaf = &pntN->grupoHijos[posicionFinal];
                    insertNode(record, pntChooseLeaf);

                    return  true;
                }
            }
        }

        // ----------- DELETE FUNCTION ----------------------

        // ----------- SEARCH FUNCTION ----------------------

        // ----------- PRINT VTK FUNCTION ----------------------
        void getAllSizesNode(Node *pntN){
            if(pntN->isLeaf == true){
                vector<float> allSizesLeaf;

                // --> Guardar bottomLeft del nodo hoja
                for(int i = 0; i < ndim; i++){
                    cout << pntN->bottomLeft[i] << " - ";
                    allSizesLeaf.push_back(pntN->bottomLeft[i]);
                }

                // --> Guardar distancias del nodo hoja (usando el upperRight)
                for(int i = 0; i < ndim; i++){
                    float aux = pntN->upperRight[i] - pntN->bottomLeft[i];
                    cout << aux  << " - " ;
                    allSizesLeaf.push_back(aux);
                }
                cout << endl;

                // --> Guardamos BottomLeft y distancias en la variable global para su impresion
                allSizes.push_back(allSizesLeaf);
                return;
            }
            else{
                cout << endl << endl << "Visit "<< pntN->grupoHijos.size() << " childrens " << endl;
                vector<float> allSizesRoot;

                // -------------- SACAMOS BL Y UR -----------------
                // --> Guardar bottomLeft del nodo padre
                for(int i = 0; i < ndim; i++){
                    cout << pntN->bottomLeft[i] << " - ";
                    allSizesRoot.push_back(pntN->bottomLeft[i]);
                }
                // --> Guardar distancias del nodo padre (usando el upperRight)
                for(int i = 0; i < ndim; i++){
                    float aux = pntN->upperRight[i] - pntN->bottomLeft[i];
                    cout << aux  << " - " ;
                    allSizesRoot.push_back(aux);
                }
                cout << endl;

                // --> Guardamos BottomLeft y distancias en la variable global para su impresion
                allSizes.push_back(allSizesRoot);

                // --> Llamamos recursivamente a la funcion con cada hijo del nodo padre
                for(int i = 0; i < pntN->grupoHijos.size(); i++){
                    cout << "  -> See child " << i+1 << ": ";
                    Node<T, ndim> *pntauxPrint;
                    pntauxPrint = &pntN->grupoHijos[i];
                    getAllSizesNode(pntauxPrint);
                }
                cout << endl;
            }
        }
};

// ################## CLASE RTREE ##############################
template< typename T, int ndim >
class RTree {
    private:
        Node<T,ndim> root;
        Node<T, ndim> *pntN;

    public:
        RTree(int _maxHijos){
            Node<T,ndim> aux(_maxHijos);
            root = aux;

            pntN = &root;
        }

        void insert(Record<T,ndim> &rec ){
            // INSERTAR
            root.insertNode(rec, pntN);
        }

        void getAllSizes(){
            root.getAllSizesNode(pntN);
        }
};