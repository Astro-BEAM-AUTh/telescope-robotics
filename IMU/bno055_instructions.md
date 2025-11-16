<H2>Βήματα για την χρήση του αισθητήρα:</H2>


<br>
Οδηγός: <a href="https://learn.adafruit.com/adafruit-bno055-absolute-orientation-sensor/">
https://learn.adafruit.com/adafruit-bno055-absolute-orientation-sensor/
</a>
<br><br><br>
1) Εγκατάσταση πακέτων κλπ με βάση τον οδηγό:
<br>
<a href="https://learn.adafruit.com/circuitpython-on-raspberrypi-linux/installing-circuitpython-on-raspberry-pi">
https://learn.adafruit.com/circuitpython-on-raspberrypi-linux/installing-circuitpython-on-raspberry-pi</a>

<br><br>
2) Ενεργοποίηση του περιβάλλοντος με την εντολή:

source env/bin/activate

<br><br>
3) Εγκατάσταση του pip3 στο περιβάλλον με τις εντολές:

python3 -m pip install --upgrade pip

python3 -m pip --version

python -m ensurepip --upgrade

(η τελευταία εντολή μάλλον δεν χρειάζεται)

<br><br>
4) Πλήρης εγκατάσταση python3:

sudo apt install python3-full

(μάλλον δεν χρειάζεται)

<br><br>
5) Εγκατάσταση της βιβλιοθήκης BNO055 με την εντολή:

pip install adafruit-circuitpython-bno055
(και όχι με την εντολή: sudo pip3 install adafruit-circuitpython-bno055)

<br><br>
6) I2C Clock Stretching
<br>
<a href="https://learn.adafruit.com/circuitpython-on-raspberrypi-linux/i2c-clock-stretching">
https://learn.adafruit.com/circuitpython-on-raspberrypi-linux/i2c-clock-stretching
</a>

<br><br>
7) Δοκιμή του προγράμματος που υπάρχει στο github (εδώ δηλαδή)

Μας ενδιαφέρει η γωνία Euler που έχει την παρακάτω μορφή:<br>
(heading, roll, pitch)

<H2>Ο αισθητήρας</H2>

BNO055 (DFRobot BNO055 9-axis Absolute Orientation Sensor)
Link to product (από το Grobotronics):
https://grobotronics.com/fermion-bno055-intelligent-9-axis-sensor.html

Περισσότερες πληροφορίες εδώ:
https://wiki.dfrobot.com/BNO055_Intelligent_9_Axis_Sensor_Module_SKU_SEN0374

Και εδώ: https://www.adafruit.com/product/2472


<H2>Χρήση</H2>

Ο αισθητήρας θα μας παρέχει το pitch και το yaw.

(https://simple.wikipedia.org/wiki/Pitch,_yaw,_and_roll)

Το yaw ουσιαστικά θα είναι η γωνία σε σχέση με τον μαγνητικό βορρά. Αν
προσθέσουμε την μαγνητική απόκλιση θα βρούμε το αζιμούθιο, δηλαδή την
γωνία σε σχέση με τον πραγματικό βορρά (ενδεχομένως φυσικά με κάποια
διαφορά).

Μαγνητική απόκλιση στην Πολυτεχνική (περίπου): 5⁰ 19’ =5,31⁰

Σχετική ίστοσελίδα: https://www.magnetic-declination.com/Greece/Athens/955650.html#google_vignette

Το pitch θα είναι η γωνία σε σχέση με το οριζόντιο επίπεδο, δηλαδή το altitude.

Στην αρχή, πριν την παρατήρηση θα παίρνουμε αυτά τα δεδομένα και θα
ξέρουμε, ουσιαστικά, που «κοιτάει» το ραδιοτηλεσκόπιο εκείνη την στιγμή, γιατί
μπορεί από την προηγουμένη παρατήρηση να έχει μετακινηθεί επίτηδες από
άνθρωπο ή από άλλη αιτία, για παράδειγμα από τον αέρα. Οπότε, θα
χρησιμοποιούμε αυτές τις συντεταγμένες ως αρχή και θα μπορούμε να
στρέψουμε το ραδιοτηλεσκόπιο προς κάποια συγκεκριμένη θέση στον ουρανό.
