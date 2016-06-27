# K2P2

Computes representation matrices for Lie groups.

# How do I get set up? #

## Clone the repo

```
git clone https://github.com/gutfeeling/K2P2.git
cd K2P2
```


## Install dependencies

1. Create a virtualenv
  - If you want to use python 2

    ```  
    virtualenv venv
    ```
  - If you want to use python 3

    ```
    virtualenv -p python3 venv
    ```

2. Activate the virtualenv

  ```
  source venv/bin/activate
  ```

3. Install the python modules
  - If you are using python 2

    ```
    pip install -r requirements.txt
    ````
  - If you are using python 3

    ```
    pip install -r requirements3.txt
    ````

## You are good to go

1. Fire up an interpreter and import the important classes

  ```python
  >>> from groups import LieGroup
  >>> from representations import Representation
  ```
