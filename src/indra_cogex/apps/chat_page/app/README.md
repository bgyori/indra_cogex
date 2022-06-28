# app

## Project setup
```sh
npm install
```

### Compiles and hot-reloads for development
```sh
npm run serve
```

### Compiles and minifies for production
```sh
npm run build
```

#### Build index.html and dist directory separately

To build public/index.html on its own, run:
```sh
npm run build-index-html
```
This will run the python script `build_index_html.py` in this directory.

To build the dist directory on its own, run:
```sh
npm run build-dist
```

### Lints and fixes files
```sh
npm run lint
```

### Customize configuration
See [Configuration Reference](https://cli.vuejs.org/config/).

## Setting up the chat app for usage with imported components from indralab-vue

### Packaging `indralab-vue` for local install
In the root of the `indralab-vue` repository:

- Run `npm run build`; this will create or update the `dist` directory (`dist` is the default name, another name might 
  be set)
- Run `npm pack`; this will package the project in tgz file that can be used as a standalone installation

### Installing indralab-vue locally

In the root directory of the Vue chat app of the `indra_cogex` repository, i.e. `indra_cogex/apps/chat_page/app`,
install the `indralab-vue` from the .tgz file. _Note that the version number might be different._

```sh
npm install /path/to/indralab-vue/indralab-vue-0.0.1.tgz
```

### Importing components from `indralab-vue`

It seems that when components from the library are nested, the nested/child components are not resolved. It's unclear 
why, but a workaround is to register all the necessary components globally in the app setup in `main.js`. When 
components are registered globally, there is no need to import them in the component files.

See more at https://vuejs.org/api/application.html#app-component

```js
// indra_cogex/src/indra_cogex/apps/chat_page/app/src/main.js
import { createApp } from "vue";

import { Statement } from "indralab-vue";
import { Evidence } from "indralab-vue";
import { CurationRow } from "indralab-vue";
import { RefLink } from "indralab-vue";
import { SourceDisplay } from "indralab-vue";
import App from "./App.vue";

// ... other setup options

const app = createApp(App);

app.component("Statement", Statement);
app.component("Evidence", Evidence);
app.component("RefLink", RefLink);
app.component("SourceDisplay", SourceDisplay);
app.component("CurationRow", CurationRow);


app.mount("#app");
```
